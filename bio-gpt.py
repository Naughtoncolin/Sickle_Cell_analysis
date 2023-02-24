# Implementation of bio-gpt

import torch
from transformers import BioGptTokenizer, BioGptForCausalLM, set_seed

tokenizer = BioGptTokenizer.from_pretrained("microsoft/biogpt")
model = BioGptForCausalLM.from_pretrained("microsoft/biogpt")

sentence = "sickle eQTL cell cycle because"
inputs = tokenizer(sentence, return_tensors="pt")

set_seed(42)

with torch.no_grad():
    beam_output = model.generate(**inputs,
                                min_length=100,
                                max_length=1024,
                                num_beams=5,
                                early_stopping=True
                                )
tokenizer.decode(beam_output[0], skip_special_tokens=True)
