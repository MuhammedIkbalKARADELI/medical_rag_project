# BU DOSYADA: LLM ile cevap üretme
from transformers import AutoTokenizer, AutoModelForCausalLM

llm_tokenizer = AutoTokenizer.from_pretrained("microsoft/BioGPT")  
llm_model = AutoModelForCausalLM.from_pretrained("microsoft/BioGPT")

def generate_answer(prompt: str) -> str:
    inputs = llm_tokenizer(prompt, return_tensors="pt")
    outputs = llm_model.generate(
        **inputs,
        max_new_tokens=300,
        do_sample=True,
        top_k = 13,
        top_p=0.65,
        temperature=0.8
    )
    result = llm_tokenizer.decode(outputs[0], skip_special_tokens=True)
    # Prompt sonrası "Answer:" kısmından kırpabilirsin:
    return result.split("Answer:")[-1].strip()