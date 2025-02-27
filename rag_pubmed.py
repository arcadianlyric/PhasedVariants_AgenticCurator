
"""
RAG with PubMed and bioRxiv
"""

from Bio import Entrez
import requests
from bs4 import BeautifulSoup
import time

# pubmed API Entrez
Entrez.email = "yingcai1123@gmail.com"  # 替换为你的邮箱

def search_pubmed(query, max_results=5):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def get_fulltext_links(pmids):
    links = []
    for pmid in pmids:
        handle = Entrez.elink(dbfrom="pubmed", id=pmid, linkname="pubmed_pmc")
        record = Entrez.read(handle)
        handle.close()
        for linkset in record[0].get("LinkSetDb", []):
            for link in linkset.get("Link", []):
                pmc_id = link["Id"]
                links.append(f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_id}/")
    return links

def fetch_abstracts(pmids):
    handle = Entrez.efetch(db="pubmed", id=pmids, rettype="abstract", retmode="text")
    abstracts = handle.read()
    handle.close()
    return abstracts

# 查询 BRCA1 和变体
# query = "BRCA1 43000604[All Fields] OR BRCA1 frameshift variant cancer"
# pmids = search_pubmed(query)

## crucial variants will be in abstracts
# abstracts = fetch_abstracts(pmids)
# print("Retrieved Abstracts:")
# print(abstracts)

## 通过 Entrez.elink 获取 PubMed 文章的 DOI 或 URL，然后尝试访问全文
# links = get_fulltext_links(pmids)
# print("PubMed Full-Text Links:")
# for link in links:
#     print(link)

## bioRxiv
"""
bioRxiv 不提供直接的 API，但可以通过网页搜索或自定义爬取全文。
使用 Python 的 requests 和 BeautifulSoup 从 bioRxiv 检索并解析预印本内容
"""

import requests
from bs4 import BeautifulSoup

def search_biorxiv(query):
    # 构造搜索 URL
    url = f"https://www.biorxiv.org/search/{query.replace(' ', '%20')}"
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
    }
    
    # 发送请求
    try:
        response = requests.get(url, headers=headers, timeout=10)
        print(f"Status Code: {response.status_code}")
        if response.status_code != 200:
            print(f"Failed to access bioRxiv: {response.status_code}")
            return []
    except requests.exceptions.RequestException as e:
        print(f"Request error: {e}")
        return []
    
    # 解析页面
    soup = BeautifulSoup(response.text, "html.parser")
    
    # 检查原始 HTML（调试用）
    with open("biorxiv_response.html", "w") as f:
        f.write(response.text)  # 保存响应以手动检查
    
    # 查找文章标题链接
    articles = soup.find_all("a", class_="highwire-cite-title")
    if not articles:
        print("No articles found. Checking alternative selectors...")
        # 尝试更通用选择器
        articles = soup.select("a[href^='/content/']")  # 匹配 /content/ 开头的链接
    
    links = ["https://www.biorxiv.org" + a["href"] for a in articles[:3] if a.get("href")]
    if not links:
        print(f"No matching links found for query: {query}")
    
    return links


def fetch_biorxiv_fulltext(url):
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
    }
    
    # 发送请求
    try:
        response = requests.get(url, headers=headers, timeout=10)
        print(f"Status Code for {url}: {response.status_code}")
        if response.status_code != 200:
            return f"Failed to fetch full text: HTTP {response.status_code}"
    except requests.exceptions.RequestException as e:
        return f"Request error: {e}"
    
    # 解析 HTML
    soup = BeautifulSoup(response.text, "html.parser")
    
    # 保存原始 HTML 以供检查
    with open("biorxiv_page.html", "w") as f:
        f.write(response.text)
    print("HTML saved to biorxiv_page.html for inspection")
    
    # 尝试原始选择器
    content = soup.find("div", class_="article-full-text")
    if content:
        return content.get_text(strip=True)
    
    # 如果未找到，尝试其他选择器
    print("Original selector 'article-full-text' not found. Trying alternatives...")
    # 常见全文容器
    alternatives = [
        soup.find("div", class_="fulltext-view"),  # 可能的全文类名
        soup.find("section", class_="article-body"),  # 文章主体
        soup.find("div", id="content-block"),  # ID 选择器
        soup.find("div", class_="main-content")  # 通用内容容器
    ]
    
    for alt in alternatives:
        if alt:
            print(f"Found content with selector: {alt.name}[class={alt.get('class')}]")
            return alt.get_text(strip=True)
    
    # 如果仍未找到，提取所有正文
    paragraphs = soup.find_all("p")
    if paragraphs:
        print("Falling back to all <p> tags")
        return "\n".join(p.get_text(strip=True) for p in paragraphs)
    
    return "No full text found after all attempts."

# test fetch_biorxiv_fulltext
# url = "https://www.biorxiv.org/content/10.1101/2021.01.19.427318v1"  
# fulltext = fetch_biorxiv_fulltext(url)
# print("Full Text Excerpt:")
# print(fulltext[:500] if fulltext else fulltext)

## test search_biorxiv
variant = "ENSG00000012048"
query = f"{variant} site:biorxiv.org"
biorxiv_links = search_biorxiv(query)
print("Search Results:")
for link in biorxiv_links:
    # print(link)
    print(f"Found {variant} in {link}")
    fulltext = fetch_biorxiv_fulltext(link)
    print(fulltext[:500])
    time.sleep(5)

