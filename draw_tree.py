
import pandas as pd


df = pd.read_csv("../experiment2/2.csv") 

len(df.origins), len(set(df.origins))

from collections import defaultdict

d = defaultdict(list)

for parent,child in zip(df.origins, df.uids):
    d[parent].append(child)

birthdays = dict(zip(df.uids, df.birthdays))

def appbirthday(x):
    return birthdays.get(x,0)

df['deliveryage'] = df.birthdays - df.origins.map(appbirthday)
deliveryage = dict(zip(df.uids, df.deliveryage))

def delivery_age(child):
    child_birthday = df[df.uids == child].birthdays
    parent = df[df.uids == child].origins
    parent_birthday = df[df.uids == parent].birthdays
    return str(child_birthday - parent_birthday)


initial = list(df[df.origins == -1].uids)

total_children = {}

def fill_total_children(parent):
    if parent not in total_children:
        total_children[parent] = 1
        
    for child in d[parent]:
        if child not in total_children:
            fill_total_children(child)
        total_children[parent] += total_children[child]

for x in initial:
    fill_total_children(x)
    
# print(total_children)

def get_newick(l):
    return "("+",".join(get_newick(d[x]) for x in l if len(d[x])>0)+")"
    
def get_newick_with_d(l):
    newl = [x for x in l if total_children.get(x,0)>100]
    if len(newl) == 0:
        return ""
    return "("+",".join(get_newick_with_d(d[x])+str(x)+":"+str(deliveryage[x]) for x in newl)+")"

import sys
sys.setrecursionlimit(1500)

newick = get_newick_with_d(initial)

with open("draw_tree.newick","w") as f:
    f.write(newick+";")


from ete3 import Tree, TreeStyle

t = Tree( "draw_tree.newick", format=1 )
circular_style = TreeStyle()
circular_style.mode = "c" # draw tree in circular mode
circular_style.scale = 20
circular_style.show_leaf_name = False

t.render("draw_tree.png", tree_style=circular_style)