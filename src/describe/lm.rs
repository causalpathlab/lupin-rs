//! The small local language model `lupin describe --lm` writes with: Qwen3
//! (Qwen Team, *Qwen3 Technical Report*, arXiv:2505.09388, 2025), run
//! through candle from the Hugging Face Hub at a pinned revision, or from a
//! local directory. Nothing is sent anywhere: the Hub is asked for the files
//! only when they are not cached yet.
//!
//! The model does two things, behind [`Lm`]: it scores a closed set of
//! continuations of a prompt, and it continues a prompt greedily until a
//! newline. Both share the prompt's KV cache, so the label choices and the
//! rationale after them pay for the prompt once.

use crate::gene_text::cli::ComputeDevice;
use anyhow::{Context, Result};
use candle_core::{DType, Device, Tensor, D};
use candle_nn::{Linear, Module, VarBuilder};
use candle_transformers::models::qwen3;
use log::info;
use std::path::{Path, PathBuf};
use tokenizers::Tokenizer;

/// What a model does for `describe`; [`Qwen3`] is the real one, the tests
/// use their own.
pub trait Lm {
    /// The mean log-probability per token of each continuation after
    /// `prompt`.
    fn score(&mut self, prompt: &str, continuations: &[String]) -> Result<Vec<f32>>;

    /// `prompt` continued greedily up to a newline or `max_tokens`, never
    /// with a token `bans` rules out; the text before the newline, its end
    /// trimmed.
    fn generate(&mut self, prompt: &str, max_tokens: usize, bans: &Bans) -> Result<String>;
}

/// The signs a sentence in words has no use for: formulas and markup.
pub const MARKUP: &str = "$\\~=*^{}<>|#_";

/// What a generation may not write.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Bans {
    /// The first token of each (with or without a space before it).
    pub words: Vec<String>,
    /// Every token with a digit or a [`MARKUP`] sign, but those of `keep`.
    pub numbers: bool,
    /// Names whose tokens stay allowed (`CD8` needs its digit).
    pub keep: Vec<String>,
}

/// Said when the model is neither cached nor at hand and the Hub cannot be
/// reached.
pub const OFFLINE: &str = "consensus needs the model: run once online or give --lm-dir";

/// The models `--lm-model` offers, each at a pinned revision (all Apache-2.0).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum LmModel {
    #[default]
    #[value(name = "qwen3-1.7b")]
    Qwen3_1_7b,
    /// A third the size, the same architecture and tokenizer: for a CPU.
    #[value(name = "qwen3-0.6b")]
    Qwen3_0_6b,
}

impl LmModel {
    #[must_use]
    pub fn id(self) -> &'static str {
        match self {
            Self::Qwen3_1_7b => "Qwen/Qwen3-1.7B",
            Self::Qwen3_0_6b => "Qwen/Qwen3-0.6B",
        }
    }

    #[must_use]
    pub fn revision(self) -> &'static str {
        match self {
            Self::Qwen3_1_7b => "70d244cc86ccca08cf5af4e1e306ecf908b1ad5e",
            Self::Qwen3_0_6b => "c1899de289a04d12100db370d81485cdf75e47ca",
        }
    }

    /// Its `--lm-model` name.
    #[must_use]
    pub fn arg(self) -> &'static str {
        match self {
            Self::Qwen3_1_7b => "qwen3-1.7b",
            Self::Qwen3_0_6b => "qwen3-0.6b",
        }
    }
}

/// The device `--device` defaults to: the GPU in a CUDA build.
#[must_use]
pub fn default_device() -> ComputeDevice {
    if cfg!(feature = "cuda") {
        ComputeDevice::Cuda
    } else {
        ComputeDevice::Cpu
    }
}

/// Which model, from where, on what.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LmSpec {
    pub model: LmModel,
    /// A directory with `config.json`, `tokenizer.json` and the
    /// `*.safetensors`; skips the Hub.
    pub local_dir: Option<PathBuf>,
    pub device: ComputeDevice,
}

impl LmSpec {
    /// The model's name as records carry it.
    #[must_use]
    pub fn name(&self) -> String {
        match &self.local_dir {
            Some(d) => format!("local:{}", d.display()),
            None => self.model.id().to_string(),
        }
    }

    /// The revision records carry: the pinned commit, or `local`.
    #[must_use]
    pub fn revision(&self) -> String {
        match &self.local_dir {
            Some(_) => "local".into(),
            None => self.model.revision().to_string(),
        }
    }
}

struct ModelFiles {
    config: PathBuf,
    tokenizer: PathBuf,
    weights: Vec<PathBuf>,
}

/// The weight files a repository splits its model into, by its
/// `model.safetensors.index.json`; one file without an index.
fn shards(index: Option<&Path>) -> Result<Vec<String>> {
    let Some(index) = index else {
        return Ok(vec!["model.safetensors".into()]);
    };
    let text =
        std::fs::read_to_string(index).with_context(|| format!("reading {}", index.display()))?;
    let v: serde_json::Value = serde_json::from_str(&text)?;
    let mut files: Vec<String> = v["weight_map"]
        .as_object()
        .context("the weight index has no weight_map")?
        .values()
        .filter_map(|f| f.as_str().map(str::to_string))
        .collect();
    files.sort();
    files.dedup();
    Ok(files)
}

/// The model's files: from the local directory, else the Hub's cache
/// (no network), else the Hub unless offline ([`crate::hub`]).
fn resolve_files(spec: &LmSpec) -> Result<ModelFiles> {
    if let Some(dir) = &spec.local_dir {
        let f = |name: &str| -> Result<PathBuf> {
            let p = dir.join(name);
            anyhow::ensure!(p.exists(), "{}: missing {name}", dir.display());
            Ok(p)
        };
        let mut weights: Vec<PathBuf> = std::fs::read_dir(dir)
            .with_context(|| format!("reading {}", dir.display()))?
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter(|p| p.extension().is_some_and(|x| x == "safetensors"))
            .collect();
        weights.sort();
        anyhow::ensure!(!weights.is_empty(), "{}: no *.safetensors", dir.display());
        return Ok(ModelFiles {
            config: f("config.json")?,
            tokenizer: f("tokenizer.json")?,
            weights,
        });
    }
    let model = crate::hub::Model::new(spec.model.id(), spec.model.revision());
    let cached = (|| {
        let index = model.cached("model.safetensors.index.json");
        let weights = shards(index.as_deref())
            .ok()?
            .iter()
            .map(|f| model.cached(f))
            .collect::<Option<Vec<_>>>()?;
        Some(ModelFiles {
            config: model.cached("config.json")?,
            tokenizer: model.cached("tokenizer.json")?,
            weights,
        })
    })();
    if let Some(files) = cached {
        return Ok(files);
    }
    let fetch = || -> Result<ModelFiles> {
        let hub = model.online(!crate::progress::enabled())?;
        let index = hub.get("model.safetensors.index.json").ok();
        let weights = shards(index.as_deref())?
            .iter()
            .map(|f| hub.get(f))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(ModelFiles {
            config: hub.get("config.json")?,
            tokenizer: hub.get("tokenizer.json")?,
            weights,
        })
    };
    fetch().map_err(|e| anyhow::anyhow!("{OFFLINE} ({e:#})"))
}

/// Qwen3 for causal decoding: the base model (whose KV cache a prompt
/// leaves behind is kept and cloned per continuation) and its head.
pub struct Qwen3 {
    /// With an empty cache.
    model: qwen3::Model,
    head: Linear,
    tokenizer: Tokenizer,
    device: Device,
    /// The last prompt run: its tokens, the model after them and the next
    /// token's log-probabilities.
    last: Option<(Vec<u32>, qwen3::Model, Vec<f32>)>,
    stops: Vec<u32>,
    /// [`Self::numeric_tokens`], once found.
    numeric: Option<Vec<u32>>,
}

impl Qwen3 {
    pub fn load(spec: &LmSpec) -> Result<Self> {
        let files = resolve_files(spec)?;
        let device = match spec.device {
            ComputeDevice::Cpu => Device::Cpu,
            ComputeDevice::Cuda => Device::new_cuda(0)?,
            ComputeDevice::Metal => Device::new_metal(0)?,
        };
        // Half precision where the GPU has it; f32 on the CPU.
        let dtype = if device.is_cuda() {
            DType::BF16
        } else {
            DType::F32
        };
        let config: qwen3::Config = serde_json::from_str(
            &std::fs::read_to_string(&files.config)
                .with_context(|| format!("reading {}", files.config.display()))?,
        )
        .context("parsing the model config")?;
        let tokenizer = Tokenizer::from_file(&files.tokenizer)
            .map_err(|e| anyhow::anyhow!("loading tokenizer: {e}"))?;
        // SAFETY: the safetensors files are memory-mapped read-only for the
        // lifetime of the loaded weights; they are not modified while mapped.
        let vb = unsafe { VarBuilder::from_mmaped_safetensors(&files.weights, dtype, &device)? };
        let model = qwen3::Model::new(&config, vb.clone()).context("loading the Qwen3 weights")?;
        let shape = (config.vocab_size, config.hidden_size);
        let head = if config.tie_word_embeddings {
            vb.pp("model.embed_tokens").get(shape, "weight")?
        } else {
            vb.pp("lm_head").get(shape, "weight")?
        };
        let stops = ["<|im_end|>", "<|endoftext|>"]
            .iter()
            .filter_map(|t| tokenizer.token_to_id(t))
            .collect();
        info!(
            "describe: {}@{} on {:?}",
            spec.name(),
            spec.revision(),
            device
        );
        Ok(Self {
            model,
            head: Linear::new(head, None),
            tokenizer,
            device,
            last: None,
            stops,
            numeric: None,
        })
    }

    /// Every token whose text has a digit or a [`MARKUP`] sign, found once.
    fn numeric_tokens(&mut self) -> Vec<u32> {
        if let Some(t) = &self.numeric {
            return t.clone();
        }
        let n = self.tokenizer.get_vocab_size(true) as u32;
        let t: Vec<u32> = (0..n)
            .filter(|&id| {
                let text = self.tokenizer.decode(&[id], false).unwrap_or_default();
                text.chars().any(|c| c.is_numeric() || MARKUP.contains(c))
            })
            .collect();
        self.numeric = Some(t.clone());
        t
    }

    fn encode(&self, text: &str) -> Result<Vec<u32>> {
        Ok(self
            .tokenizer
            .encode(text, false)
            .map_err(|e| anyhow::anyhow!("tokenizing: {e}"))?
            .get_ids()
            .to_vec())
    }

    /// Log-probabilities over the vocabulary from `hidden` ([1, n, h]) at
    /// every position, as rows.
    fn log_probs(&self, hidden: &Tensor) -> Result<Vec<Vec<f32>>> {
        let logits = self.head.forward(hidden)?.to_dtype(DType::F32)?;
        let lp = candle_nn::ops::log_softmax(&logits, D::Minus1)?;
        Ok(lp.squeeze(0)?.to_vec2::<f32>()?)
    }

    /// `tokens` run through `model` from `offset`: the next token's
    /// log-probabilities.
    fn feed(&self, model: &mut qwen3::Model, tokens: &[u32], offset: usize) -> Result<Vec<f32>> {
        let input = Tensor::new(tokens, &self.device)?.unsqueeze(0)?;
        let hidden = model.forward(&input, offset)?;
        let last = hidden.narrow(1, tokens.len() - 1, 1)?;
        Ok(self.log_probs(&last)?.remove(0))
    }

    /// The model after `tokens`, and the next token's log-probabilities:
    /// from the last prompt when it starts `tokens`.
    fn prefill(&mut self, tokens: &[u32]) -> Result<(qwen3::Model, Vec<f32>)> {
        anyhow::ensure!(!tokens.is_empty(), "an empty prompt");
        let reuse = self
            .last
            .as_ref()
            .filter(|(t, ..)| tokens.starts_with(t))
            .map(|(t, m, lp)| (t.len(), m.clone(), lp.clone()));
        let (at, mut model, mut lp) = match reuse {
            Some(r) => r,
            None => (0, self.model.clone(), Vec::new()),
        };
        if at < tokens.len() {
            lp = self.feed(&mut model, &tokens[at..], at)?;
        }
        self.last = Some((tokens.to_vec(), model.clone(), lp.clone()));
        Ok((model, lp))
    }
}

/// The index of the largest of `v`, the first on a tie.
fn argmax(v: &[f32]) -> usize {
    let mut best = 0;
    for (i, &x) in v.iter().enumerate() {
        if x > v[best] {
            best = i;
        }
    }
    best
}

impl Lm for Qwen3 {
    fn score(&mut self, prompt: &str, continuations: &[String]) -> Result<Vec<f32>> {
        let tokens = self.encode(prompt)?;
        let (model, first) = self.prefill(&tokens)?;
        continuations
            .iter()
            .map(|c| {
                let ids = self.encode(c)?;
                anyhow::ensure!(!ids.is_empty(), "an empty continuation");
                let mut total = first[ids[0] as usize];
                if ids.len() > 1 {
                    let mut m = model.clone();
                    let input = Tensor::new(&ids[..ids.len() - 1], &self.device)?.unsqueeze(0)?;
                    let hidden = m.forward(&input, tokens.len())?;
                    for (row, &next) in self.log_probs(&hidden)?.iter().zip(&ids[1..]) {
                        total += row[next as usize];
                    }
                }
                Ok(total / ids.len() as f32)
            })
            .collect()
    }

    fn generate(&mut self, prompt: &str, max_tokens: usize, bans: &Bans) -> Result<String> {
        let tokens = self.encode(prompt)?;
        let (mut model, mut lp) = self.prefill(&tokens)?;
        let mut ban: Vec<u32> = Vec::new();
        if bans.numbers {
            let keep: Vec<u32> = (bans.keep.iter())
                .flat_map(|k| [k.clone(), format!(" {k}")])
                .map(|k| self.encode(&k))
                .collect::<Result<Vec<_>>>()?
                .concat();
            let numeric = self.numeric_tokens();
            ban.extend(numeric.iter().filter(|t| !keep.contains(t)));
        }
        for w in &bans.words {
            for form in [w.clone(), format!(" {w}")] {
                ban.extend(self.encode(&form)?.first().copied());
            }
        }
        let mut out: Vec<u32> = Vec::new();
        for step in 0..max_tokens {
            for &b in &ban {
                lp[b as usize] = f32::NEG_INFINITY;
            }
            let next = argmax(&lp) as u32;
            if self.stops.contains(&next) {
                break;
            }
            out.push(next);
            let text = self.tokenizer.decode(&out, true).unwrap_or_default();
            if text.contains('\n') || step + 1 == max_tokens {
                break;
            }
            lp = self.feed(&mut model, &[next], tokens.len() + step)?;
        }
        let text = self
            .tokenizer
            .decode(&out, true)
            .map_err(|e| anyhow::anyhow!("decoding: {e}"))?;
        // A cut through a character's bytes leaves a replacement sign.
        let line = text.split('\n').next().unwrap_or("");
        Ok(line.trim_end_matches(['\u{FFFD}', ' ']).to_string())
    }
}

#[cfg(test)]
#[path = "tests/lm.rs"]
mod tests;
