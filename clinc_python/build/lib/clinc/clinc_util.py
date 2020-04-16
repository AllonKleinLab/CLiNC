import numpy as np, os, matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
from fastcluster import linkage

def get_hierch_order(hm, dist_metric='euclidean', linkage_method='ward'):
    np.random.seed(0)
    D = pdist(hm, dist_metric)
    Z = linkage(D, linkage_method)
    n = len(Z) + 1
    cache = dict()
    for k in range(len(Z)):
        c1, c2 = int(Z[k][0]), int(Z[k][1])
        c1 = [c1] if c1 < n else cache.pop(c1)
        c2 = [c2] if c2 < n else cache.pop(c2)
        cache[n+k] = c1 + c2
    o = np.array(cache[2*len(Z)])
    return o

def load_data(input_data_path):
    barcode_counts = np.loadtxt(input_data_path, delimiter='\t', skiprows=1)
    celltype_names = open(input_data_path).readline().replace('\n','').replace('\r','').strip('\t').split('\t')
    print('Loaded',barcode_counts.shape[0],'barcodes in',barcode_counts.shape[1],'cell types')
    return barcode_counts, celltype_names
    
def plot_barcode_counts(output_directory, barcode_counts, celltype_names):
    o = get_hierch_order(barcode_counts)
    plt.imshow(np.log(barcode_counts[o,:]+1)/np.log(10), aspect='auto',cmap=plt.cm.Reds, vmax=1)
    plt.xticks(np.arange(barcode_counts.shape[1])+.4, celltype_names, rotation=70, ha='right', fontsize=12)
    plt.yticks([])
    cbar = plt.colorbar()
    cbar.set_label('Number of barcodes (log10)', rotation=270, fontsize=12, labelpad=20)
    plt.gcf().set_size_inches((4,6))
    plt.tight_layout()
    plt.savefig(output_directory+'/barcode_counts.png', dpi=300)

def make_output_dir(output_directory):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        
def get_normalized_covariance(data):
    cc = np.cov(data.T)
    mm = np.mean(data,axis=0) + .0001
    X,Y = np.meshgrid(mm,mm)
    cc = cc / X / Y
    return cc

def get_normalized_covariance_boostrap(data, n_iters=5000):
    ccs = []
    for i in range(n_iters):
        ix = np.random.randint(0,data.shape[0],data.shape[0])
        cc = np.cov(data[ix,:].T)
        mm = np.mean(data[ix,:],axis=0)
        X,Y = np.meshgrid(mm,mm)
        cc = cc / X / Y
        ccs.append(cc)
    return np.array(ccs)

def plot_normalized_covariance(output_directory, X, celltype_names):
    vmax = (np.percentile(X-np.diag(np.diag(X)),95) + np.percentile(X-np.diag(np.diag(X)),98))/2
    plt.imshow(X, vmax=vmax)
    cbar = plt.colorbar()
    cbar.set_label('Normnalized covariance', rotation=270, fontsize=12, labelpad=20)
    plt.xticks(np.arange(X.shape[0])+.4, celltype_names, ha='right', rotation=60, fontsize=10);
    plt.yticks(np.arange(X.shape[0])+.4, celltype_names, fontsize=10);
    plt.gcf().set_size_inches((6,4))
    plt.tight_layout()
    plt.savefig(output_directory+'/normalized_covariance.png', dpi=300)
    
def build_hierarchy(barcode_counts):
    X_history = []
    merged_pairs_history = []
    node_names_history = []
    node_groups = {i:[i] for i in range(barcode_counts.shape[1])}

    parent_map = {}
    barcode_counts_tmp = np.array(barcode_counts)
    node_names = list(range(barcode_counts.shape[1]))
    next_node = barcode_counts.shape[1]

    while len(node_names) > 2: 
        node_names_history.append(node_names)
        X = get_normalized_covariance(barcode_counts_tmp)
        X_history.append(np.array(X))
        floor = X.min() - 100
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                if i >= j: X[i,j] = floor
    

        ii = np.argmax(X.max(1))
        jj = np.argmax(X.max(0))
        merged_pairs_history.append((ii,jj))
        node_groups[next_node] = node_groups[node_names[ii]]+node_groups[node_names[jj]]

        parent_map[node_names[ii]] = next_node
        parent_map[node_names[jj]] = next_node

        ix = np.min([ii,jj])
        node_names = [n for n in node_names if not n in np.array(node_names)[np.array([ii,jj])]]
        new_ix = np.array([i for i in range(barcode_counts_tmp.shape[1]) if not i in [ii,jj]])

        if len(new_ix)==0: break
        new_column = barcode_counts_tmp[:,np.array([ii,jj])].sum(1)[:,None]
        barcode_counts_tmp = barcode_counts_tmp[:,new_ix]
        barcode_counts_tmp = np.hstack((barcode_counts_tmp[:,:ix], new_column, barcode_counts_tmp[:,ix:]))
        node_names.insert(ix,next_node)
        next_node += 1


    for i in node_names:
        parent_map[i] = next_node

    return parent_map, node_groups, (X_history, merged_pairs_history, node_names_history)


def plot_neighbor_joining(output_directory, node_groups, celltype_names, X_history, merged_pairs_history, node_names_history):
    fig,axs = plt.subplots(1,len(X_history))
    for i,X in enumerate(X_history):
        vmaxx = 40
        axs[i].imshow(X,vmax=vmaxx)
        ii,jj = merged_pairs_history[i]
        axs[i].scatter([jj],[ii],s=100, marker='*', c='white')

        column_groups = [node_groups[n] for n in node_names_history[i]]
        column_labels = [' + '.join([celltype_names[n] for n in grp]) for grp in column_groups]
        axs[i].set_xticks(np.arange(X.shape[1])+.2)
        axs[i].set_xticklabels(column_labels, rotation=90, ha='right')
        axs[i].set_xlim([-.5,X.shape[1]-.5])
        axs[i].set_ylim([X.shape[1]-.5,-.5])
    fig.set_size_inches((16,4))
    plt.savefig(output_directory+'/neighbor_joint_heatmaps.pdf')

def print_hierarchy(parent_map, celltype_names):    
    child_map = {i:[] for i in set(list(parent_map.values())+list(parent_map.keys()))}
    for i,j in parent_map.items():
        child_map[j].append(i)

    leaf_names = {i:n for i,n in enumerate(celltype_names)}
    def get_newick(n):
        if n in leaf_names: return leaf_names[n]
        else: return '('+','.join([get_newick(nn) for nn in sorted(child_map[n])[::-1]])+')'
    tree_string = get_newick(np.max(list(child_map.keys())))+';'
    
    from ete3 import Tree
    t = Tree(tree_string)
    print(t)
    
    
def mrca(i,j, parent_map):
    ips = set([i])
    jps = set([j])
    while len(ips.intersection(jps))==0:
        if i in parent_map:
            ips.add(parent_map[i])
            i = parent_map[i]
        if j in parent_map:
            jps.add(parent_map[j])
            j = parent_map[j]
    return list(ips.intersection(jps))[0]

def get_tree_heights(parent_map, N):
    tree_heights = {}
    for i in range(N): tree_heights[i] = 0
    for j in range(N, len(parent_map)+1):
        prev_heights = [tree_heights[i] for i in parent_map if parent_map[i]==j]
        tree_heights[j] = np.max(prev_heights)+1
    return tree_heights
            

def tree_dist(i,j, tree_heights, parent_map):
    k = mrca(i,j, parent_map)
    if k in [i,j]: return 0
    else: return tree_heights[k]

    

def detect_symmetry_violations(barcode_counts, parent_map, symmetry_violation_FDR):
    Xbootstrap = get_normalized_covariance_boostrap(barcode_counts)
    X = get_normalized_covariance(barcode_counts)
    triples = []
    diffs = []
    diffs_upper = []
    diffs_lower = []
    val_pairs = []
    for j in range(barcode_counts.shape[1]):
        for i in range(barcode_counts.shape[1]):
            for k in range(barcode_counts.shape[1]):
                if len(set([i,j,k])) ==3:
                    p1 = mrca(i,j, parent_map)
                    p2 = mrca(i,k, parent_map)
                    if p2 == mrca(p1,p2, parent_map) and p1 != p2:
                        triples.append((i,j,k))
                        diffs.append(X[j,k]-X[i,k])
                        diffs_upper.append(np.percentile(Xbootstrap[:,j,k] - Xbootstrap[:,i,k], 100*(1-symmetry_violation_FDR)))
                        diffs_lower.append(np.percentile(Xbootstrap[:,j,k] - Xbootstrap[:,i,k], symmetry_violation_FDR*100))
                        val_pairs.append((X[i,k], X[j,k]))

    threshold = np.median(np.abs(diffs))
    violations = [triples[i] for i in np.nonzero(np.array(diffs_lower) > threshold)[0]]
    print('Detected', len(violations), 'instances of symmetry violation passing a threshold of',threshold,'with FDR', symmetry_violation_FDR)
    return violations, (diffs, diffs_lower, diffs_upper, threshold, val_pairs)
    

def plot_violations(output_directory, diffs, diffs_lower, diffs_upper, threshold, val_pairs):
    fig,axs = plt.subplots(1,2)
    
    o = np.argsort(diffs)
    o = o[np.array(diffs)[o]>0]
    o1 = o[np.array(diffs_lower)[o] < threshold]
    o2 = o[np.array(diffs_lower)[o] > threshold]

    axs[0].errorbar(np.arange(len(o1)), np.array(diffs)[o1], yerr=(np.array(diffs)[o1] - np.array(diffs_lower)[o1]), c='gray', linewidth=.5)
    axs[0].scatter(np.arange(len(o1)), np.array(diffs)[o1], c='gray',s=15, edgecolor='k', linewidth=0.2)
    axs[0].errorbar(np.arange(len(o1),len(o)), np.array(diffs)[o2], yerr=(np.array(diffs)[o2] - np.array(diffs_lower)[o2]), c='gray', zorder=1, linewidth=.5)
    axs[0].scatter(np.arange(len(o1),len(o)), np.array(diffs)[o2], c='red', zorder=2,s=15, edgecolor='k', linewidth=0.2)

    axs[0].plot([0,len(o)],[threshold,threshold],'--k')
    axs[0].set_xlabel('Putative symmetric triples')
    axs[0].set_ylabel('Diff of normed covariance')
    
    val_pairs = np.array(val_pairs)
    axs[1].scatter(val_pairs[:,0], val_pairs[:,1],s=30, edgecolor='k', linewidth=0.5, c='lightgray', zorder=2)
    ff = np.array(diffs_lower) > threshold
    axs[1].scatter(val_pairs[ff,0], val_pairs[ff,1],s=30, edgecolor='k', linewidth=0.5, c='r', zorder=3)
    axs[1].set_xlabel('Normalized covariance (j,k)')
    axs[1].set_ylabel('Normalized covariance (i,k)')
    
    fig.set_size_inches((6,2.3))   
    fig.subplots_adjust(wspace=.5)
    plt.savefig(output_directory+'/symmetry_violations.pdf')

    

def get_violations(ip,jp,N, parent_map, tree_heights):
    violations = []
    if mrca(ip,jp,parent_map) in [ip,jp]: return []
    for i in range(N):
        for j in range(N):
            for k in range(N):
                if len(set([i,j,k])) ==3:
                    p1 = mrca(i,j, parent_map)
                    p2 = mrca(i,k, parent_map)
                    if p2 == mrca(p1,p2, parent_map) and p1 != p2:
                        if 0 < tree_dist(k,ip, tree_heights, parent_map) and  tree_dist(k,ip,tree_heights, parent_map) <  tree_dist(k,jp,tree_heights, parent_map):
                            if tree_dist(i,ip,tree_heights, parent_map) == 0 and tree_dist(j,ip,tree_heights, parent_map) > 0:
                                violations.append((i,j,k))
                                #violations.append((j,i,k))

                        if tree_dist(k,ip,tree_heights, parent_map) > tree_dist(k,jp,tree_heights, parent_map):
                            if tree_dist(i,ip,tree_heights, parent_map) > 0 and tree_dist(j,ip,tree_heights, parent_map) == 0:
                                violations.append((i,j,k))
                                #violations.append((j,i,k))

                        if tree_dist(k,ip,tree_heights, parent_map) == 0:                        
                            if tree_dist(i,ip,tree_heights, parent_map) >  tree_dist(i,jp,tree_heights, parent_map) and  tree_dist(i,jp,tree_heights, parent_map) > tree_dist(j,jp,tree_heights, parent_map):
                                violations.append((i,j,k))
                                #violations.append((j,i,k))
                            
    return violations

from SetCoverPy import setcover

def detect_cross_tree_transitions(parent_map, violations, N):
    tree_heights = get_tree_heights(parent_map, N)

    violation_templates = {}
    for i in parent_map:
        for j in parent_map:
            if i != j:
                violation_templates[(i,j)] = get_violations(i,j,N,parent_map, tree_heights)
                
    initial_match_scores = []
    initial_transitions = []
    for transition,vs in violation_templates.items():
        if len(vs) > 0:
            matches = len(set(violations).intersection(vs))
            initial_match_scores.append(matches-len(vs))
            initial_transitions.append(transition)   
    
    transitions = [t for t in violation_templates.keys() if len(violation_templates[t])>0]
    templates = [violation_templates[t] for t in transitions]
    templates_union = set([])
    for template in templates: templates_union = templates_union.union(template)
    violation_order = [v for v in violations if v in templates_union]
    
    a_matrix = np.zeros((len(violation_order),len(templates)))
    for j,template in enumerate(templates):
        for i,v in enumerate(violation_order):
            a_matrix[i,j] = 1 if v in template else 0
    cost = np.array([len([t for t in template if not t in violations])/len(template) for template in templates])+1
    
    g = setcover.SetCover(a_matrix>0, cost)
    solution, time_used = g.SolveSCP()
    nz = np.nonzero(g.s)[0] 
    final_transitions = [transitions[i] for i in nz if cost[i] <= 1.5]
    
    num_violations_predicted = [len(violation_templates[t]) for t in final_transitions]
    num_violations_explained = [len(set(violation_templates[t]).intersection(violations)) for t in final_transitions]
    explained = []
    for t in final_transitions: explained += [v for v in violations if v in violation_templates[t]]
    total_explained = len(set(explained))
    return final_transitions, num_violations_predicted, num_violations_explained, total_explained

def print_cross_tree_transitions(final_transitions, num_violations_predicted, num_violations_explained, total_explained, total_violations, node_groups, celltype_names):
    print(total_explained, 'out of', total_violations, 'symmetry violations can be explained by the following cross-tree transitions:\n')
    print('Explained symmetry violations    Proportion matching    Transition')
    
    labels_ord = ['+'.join([celltype_names[i] for i in node_groups[n]]) for  n in sorted(node_groups.keys())]
    for (i,j), n_explained, n_predicted in zip(final_transitions, num_violations_explained, num_violations_predicted):
        prop_match_str = repr(round(n_explained/n_predicted,2))
        explained_str = repr(n_explained)
        transition_str = labels_ord[j]+' -> '+labels_ord[i]
        print(explained_str+' '*(33-len(explained_str))+prop_match_str+' '*(23-len(prop_match_str))+transition_str)
        

