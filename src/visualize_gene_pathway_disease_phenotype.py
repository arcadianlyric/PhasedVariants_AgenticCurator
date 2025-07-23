import json
import networkx as nx
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict

def load_json(INPUT):
    """load JSON from previous step"""
    with open(INPUT, 'r') as f:
        data = json.load(f)
    return data

def build_network_graph(data):
    """construct networkX graph"""
    G = nx.Graph()
    
    # add nodes
    genes = data.get('genes', [])
    for gene in genes:
        G.add_node(gene, type='gene', size=10)
    
    # print("Common Pathways:", data.get('common_pathways', []))
    # print("Common Diseases:", data.get('common_diseases', []))
    # print("Common Phenotypes:", data.get('common_phenotypes', []))
    # print("Details:", data.get('details', {}))
    
    # add nodes and edges
    for pathway in data.get('common_pathways', []):
        if isinstance(pathway, dict) and 'name' in pathway:
            name = pathway['name']
            associated_genes = pathway.get('associated_genes', [])
            size = 15 + 5 * len(associated_genes)
            G.add_node(name, type='pathway', size=size, gene_count=len(associated_genes))
            for gene in associated_genes:
                G.add_edge(gene, name)
        else:
            print(f"Skipping invalid pathway: {pathway}")
    
    # add disease 
    disease_genes = defaultdict(set)
    for gene in data.get('details', {}):
        diseases = data['details'][gene].get('diseases', [])
        for disease in diseases:
            if isinstance(disease, str):
                disease_genes[disease].add(gene)
                size = 15 + 5 * len(disease_genes[disease])
                G.add_node(disease, type='disease', size=size, gene_count=len(disease_genes[disease]))
                G.add_edge(gene, disease)
            else:
                print(f"Skipping invalid disease for {gene}: {disease}")
    
    # add phenotype
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
    """interactive plot with Plotly """
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
                        title='Gene-Pathway-Disease-Phenotype Network (All Diseases)',
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
                    ))
    
    fig.write_html(output_html)
    print(f"Network graph saved as {output_html}")

def plot_bar_chart(data, output_png='importance_bar.png'):
    """barplot"""
    items = []

    for p in data.get('common_pathways', []):
        if isinstance(p, dict) and 'name' in p:
            items.append({'name': p['name'], 'type': 'pathway', 'gene_count': len(p.get('associated_genes', []))})

    disease_genes = defaultdict(set)
    for gene in data.get('details', {}):
        for d in data['details'][gene].get('diseases', []):
            if isinstance(d, str):
                disease_genes[d].add(gene)
                items.append({'name': d, 'type': 'disease', 'gene_count': len(disease_genes[d])})

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
    plt.title('Gene Coverage by Pathway, Disease, and Phenotype (All Diseases)')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    plt.savefig(output_png)
    plt.close()
    print(f"Bar chart saved as {output_png}")

def plot_graph():
    data = load_json(INPUT)
    
    G = build_network_graph(data)

    plot_network_graph(G, output_html=OUTPUT)
    
    # plot_bar_chart(data, output_png='importance_bar.png')

if __name__ == "__main__":
    INPUT = f"../results/gene_associations.json"
    OUTPUT = f'../results/network_graph.html'
    plot_graph()