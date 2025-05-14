## 
'''
input:
gene name
pubmed abstract from pubmed.py
output:
Q: summarize the function of the {keyword} 
P2RX5 is a P2X receptor involved in T cell activation and immunoregulation. It also plays a role in brown adipocyte differentiation and energy homeostasis. Additionally, it functions as an ATP-gated cation channel mediating physiological processes.
Q: what is the consequence of loss of function of {keyword}
missing ATP mechanism
Loss of P2RX5 function impairs T cell activation, leading to increased IL-10 production and altered immunoregulation. It also reduces bone loss in periodontitis by downregulating pro-inflammatory cytokines like IL1b and IL17a. Additionally, P2RX5 loss hampers brown adipocyte differentiation and energy homeostasis, potentially contributing to metabolic dysregulation.
Q: what is the consequence of loss of function of {2 genes}
missing P2RX5
The context provided does not mention P2RX5, so its loss of function consequences are not addressed here. CYP2D6 loss of function (e.g., in poor metabolizers) leads to reduced metabolism of drugs like antidepressants and opioids, affecting therapeutic outcomes. Further studies are needed to clarify interactions between CYP2D6 and other genes like P2RX5.

RAG: use input.pdf to answer questions

conda create -n llm (tensorflow incompatible)
pip install langchain 
pip install langchain_community
pip install langchain_chroma
pip install langchain_ollama
pip instal pdfplumber

ollama pull deepseek-r1:7b
ollama pull qwen3:8b
ollama pull nomic-embed-text
'''

## process input.pdf
from langchain_community.document_loaders import PDFPlumberLoader, TextLoader
from langchain.text_splitter import RecursiveCharacterTextSplitter

# input_text = '/Users/yc/Documents/GitHub/llm_genetics_assistant/ENSG00000012048.pdf'
input_text ='pubmed_response.txt'
keyword = 'P2RX5 and CYP2D6'
keyword = keyword.lower()


if input_text.endswith('.pdf'):
    loader=PDFPlumberLoader(input_text)
elif input_text.endswith('.txt'):
    loader=TextLoader(input_text)
docs =loader.load()


text_splitter = RecursiveCharacterTextSplitter(chunk_size=500,chunk_overlap=0)
all_splits = text_splitter.split_documents(docs)

## text splitting 
from langchain_chroma import Chroma
from langchain_ollama import OllamaEmbeddings

local_embeddings = OllamaEmbeddings(model="nomic-embed-text")

vectorstore = Chroma.from_documents(documents=all_splits, embedding=local_embeddings)

## create embedding database
from langchain_core.output_parsers import StrOutputParser
from langchain_core.prompts import ChatPromptTemplate
from langchain_ollama import ChatOllama

model = ChatOllama(
    model="qwen3:8b",
)

# prompt = ChatPromptTemplate.from_template(
#     f"Summarize the function of {keyword} in these retrieved docs: {docs}"
# )

# create chain 
def format_docs(docs):
    return "\n\n".join(doc.page_content for doc in docs)

# chain = {"docs": format_docs} | prompt | model | StrOutputParser()
# question = "What is the purpose of the input doc ?"
# docs = vectorstore.similarity_search(question)
# chain.invoke(docs)

## retrival and answer
from langchain_core.runnables import RunnablePassthrough

RAG_TEMPLATE = """
You are an assistant for question-answering tasks. Use the following pieces of retrieved context to answer the question. If you don't know the answer, just say that you don't know. Use three sentences maximum and keep the answer concise.

<context>
{context}
</context>

Answer the following question:

{question}"""

rag_prompt = ChatPromptTemplate.from_template(RAG_TEMPLATE)
## llm 的上下文限制
retriever = vectorstore.as_retriever(search_kwargs={"k": 10})

qa_chain = (
    {"context": retriever | format_docs, "question": RunnablePassthrough()}
    | rag_prompt
    | model
    | StrOutputParser()
)

# question = f"summarize the function of {keyword} based on the provided context"
question = f"what is the consequence of loss of function of {keyword}"
answer = qa_chain.invoke(question)
print(answer)
# Loss of P2RX5 function impairs T cell activation, leading to increased IL-10 production and altered immunoregulation. It also reduces bone loss in periodontitis by downregulating pro-inflammatory cytokines like IL1b and IL17a. Additionally, P2RX5 loss hampers brown adipocyte differentiation and energy homeostasis, potentially contributing to metabolic dysregulation.