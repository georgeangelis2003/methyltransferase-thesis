from ete4 import Tree
import random

tree = Tree("/home/angelis/thesis/reference_proteomes/new_filtered_diamond_tsv/iqtree_LG_G_Tree.treefile")

def get_organism_name(node_name):
    if node_name and "_" in node_name:
        return node_name.split("_")[-1]
    return "Unknown"

all_orgs = list(set(get_organism_name(leaf.name) for leaf in tree.leaves()))
org_to_color = {org: f"#{random.randint(0, 0xFFFFFF):06x}" for org in all_orgs}

for leaf in tree.leaves():
    org = get_organism_name(leaf.name)
    color = org_to_color.get(org, "#000000")
    
    leaf.props["container_color"] = color
    leaf.props["hz_line_color"] = color
    leaf.props["vt_line_color"] = color


tree.render_sm("my_massive_tree.html")

print(f"Done! Created an interactive map for {len(all_orgs)} organisms.")
print("Download 'my_massive_tree.html' and open it in your browser.")