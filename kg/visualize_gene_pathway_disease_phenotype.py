import json
import networkx as nx
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import pandas as pd

def load_json(json_path):
    """加载JSON文件"""
    with open(json_path, 'r') as f:
        data = json.load(f)
    return data

def build_network_graph(data):
    """构建NetworkX图"""
    G = nx.Graph()
    
    # 添加基因节点
    genes = data.get('genes', [])
    for gene in genes:
        G.add_node(gene, type='gene', size=10)
    
    # 调试：打印JSON结构
    # print("Common Pathways:", data.get('common_pathways', []))
    # print("Common Diseases:", data.get('common_diseases', []))
    # print("Common Phenotypes:", data.get('common_phenotypes', []))
    
    # 添加路径节点和边
    for pathway in data.get('common_pathways', []):
        if isinstance(pathway, dict) and 'name' in pathway:
            name = pathway['name']
            associated_genes = pathway.get('associated_genes', [])
            size = 15 + 5 * len(associated_genes)  # Size based on gene count
            G.add_node(name, type='pathway', size=size, gene_count=len(associated_genes))
            for gene in associated_genes:
                G.add_edge(gene, name)
        else:
            print(f"Skipping invalid pathway: {pathway}")
    
    # 添加疾病节点和边
    for disease in data.get('common_diseases', []):
        if isinstance(disease, dict) and 'name' in disease:
            name = disease['name']
            associated_genes = disease.get('associated_genes', [])
            size = 15 + 5 * len(associated_genes)
            G.add_node(name, type='disease', size=size, gene_count=len(associated_genes))
            for gene in associated_genes:
                G.add_edge(gene, name)
        else:
            print(f"Skipping invalid disease: {disease}")
    
    # 添加表型节点和边
    for phenotype in data.get('common_phenotypes', []):
        if isinstance(phenotype, dict) and 'name' in phenotype:
            name = phenotype['name']
            associated_genes = phenotype.get('associated_genes', [])
            size = 15 + 5 * len(associated_genes)
            G.add_node(name, type='phenotype', size=size, gene_count=len(associated_genes))
            for gene in associated_genes:
                G.add_edge(gene, name)
        else:
            print(f"Skipping invalid phenotype: {phenotype}")
    
    return G

def plot_network_graph(G, output_html='network_graph.html'):
    """使用Plotly绘制交互式网络图"""
    pos = nx.spring_layout(G)
    
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
    
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines'
    )
    
    node_x = []
    node_y = []
    node_sizes = []
    node_colors = []
    node_text = []
    
    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_sizes.append(G.nodes[node]['size'])
        node_type = G.nodes[node]['type']
        color = {'gene': 'blue', 'pathway': 'green', 'disease': 'red', 'phenotype': 'purple'}[node_type]
        node_colors.append(color)
        text = f"Name: {node}<br>Type: {node_type}"
        if node_type != 'gene':
            text += f"<br>Gene Count: {G.nodes[node]['gene_count']}"
        node_text.append(text)
    
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        text=[n[:10] for n in G.nodes()],
        textposition='top center',
        hoverinfo='text',
        hovertext=node_text,
        marker=dict(
            showscale=False,
            color=node_colors,
            size=node_sizes,
            line_width=2
        )
    )
    
    fig = go.Figure(data=[edge_trace, node_trace],
                    layout=go.Layout(
                        title='Gene-Pathway-Disease-Phenotype Network',
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
                    ))
    
    fig.write_html(output_html)
    print(f"Network graph saved as {output_html}")

def plot_bar_chart(data, output_png='importance_bar.png'):
    """绘制柱状图展示基因覆盖度"""
    items = []
    for p in data.get('common_pathways', []):
        if isinstance(p, dict) and 'name' in p:
            items.append({'name': p['name'], 'type': 'pathway', 'gene_count': len(p.get('associated_genes', []))})
    for d in data.get('common_diseases', []):
        if isinstance(d, dict) and 'name' in d:
            items.append({'name': d['name'], 'type': 'disease', 'gene_count': len(d.get('associated_genes', []))})
    for ph in data.get('common_phenotypes', []):
        if isinstance(ph, dict) and 'name' in ph:
            items.append({'name': ph['name'], 'type': 'phenotype', 'gene_count': len(ph.get('associated_genes', []))})
    
    df = pd.DataFrame(items)
    if df.empty:
        print("No valid data to plot for bar chart.")
        return
    
    df = df.sort_values('gene_count', ascending=False)
    
    plt.figure(figsize=(12, 6))
    colors = df['type'].map({'pathway': 'green', 'disease': 'red', 'phenotype': 'purple'})
    plt.bar(df['name'], df['gene_count'], color=colors)
    plt.xlabel('Pathway/Disease/Phenotype')
    plt.ylabel('Gene Count')
    plt.title('Gene Coverage by Pathway, Disease, and Phenotype')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    plt.savefig(output_png)
    plt.close()
    print(f"Bar chart saved as {output_png}")

def plot_graph():
    # JSON文件路径
    json_path = "gene_associations.json"
    
    # 加载JSON数据
    data = load_json(json_path)
    
    # 构建网络图
    G = build_network_graph(data)
    
    # 绘制交互式网络图
    plot_network_graph(G, output_html='network_graph.html')
    
    # 绘制柱状图
    # plot_bar_chart(data, output_png='importance_bar.png')

if __name__ == "__main__":
    plot_graph()