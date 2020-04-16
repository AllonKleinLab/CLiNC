#' @export
make_output_dir <- function(output_directory) {
    dir.create(output_directory, showWarnings = FALSE)
}

#' @export
load_data <- function(input_data_path) {
    barcode_counts <- read.table(file = input_data_path, sep = '\t', header = TRUE)
    print(paste0('Loaded ',nrow(barcode_counts),' barcodes in ',ncol(barcode_counts),' cell types'))
    return(barcode_counts)
}

#' @export
get_normalized_covariance <- function(barcode_counts) {
    cc <- cov(barcode_counts)
    mm <- colMeans(data.matrix(barcode_counts))
    return(cc / mm[row(cc)] / mm[col(cc)])
}

#' @export
plot_normalized_covariance <- function(output_directory, X) {
    vmax <- (quantile(X - diag(diag(X)),.95) + quantile(X - diag(diag(X)),.98))/2
    pdf(paste0(output_directory,'/normalized_covariance.pdf'))
    heatmap3::heatmap3(pmin(X,vmax), Rowv=NA, Colv=NA, 
            col = colorRampPalette(c("white", "black"))(1024),
            main="Normalized covariance")  
    dev.off()
}

#' @export
build_hierarchy <- function(barcode_counts) {
    X_history <- list()
    merged_pairs_history <- list()
    node_names_history <- list()
    node_groups <- list()
    parent_map <- list()
    for (i in 1:ncol(barcode_counts)){
        node_groups[[i]] <- list(i)
    }
    
    barcode_counts_tmp <- data.matrix(barcode_counts)
    node_names <- 1:ncol(barcode_counts)
    next_node <- ncol(barcode_counts)+1
    
    while (length(node_names) > 2) {
        node_names_history <- c(node_names_history, list(node_names))
        X <- get_normalized_covariance(barcode_counts_tmp)
        X_history <- c(X_history, list(X))
        
        Xfloor = min(X)-100
        for (i in 1:nrow(X)) {
            for (j in 1:ncol(X)) {
                if (i >= j) {
                    X[i,j] <- Xfloor
                }
            }
        }

        ii <- ramify::argmax(t(data.matrix(matrixStats::rowMaxs(X))))
        jj <- ramify::argmax(t(data.matrix(matrixStats::colMaxs(X))))
        merged_pairs_history <- c(merged_pairs_history,list(c(ii,jj)))
        node_groups[[next_node]] <- c(node_groups[[node_names[[ii]]]] , node_groups[[node_names[[jj]]]])
        parent_map[[node_names[[ii]]]] <- next_node
        parent_map[[node_names[[jj]]]] <- next_node
        
        ix <- min(ii,jj)
        node_names <- unlist(comprehenr::to_list(for(n in node_names) if(n != node_names[[ii]] & n != node_names[[jj]]) n))
        new_ix <- comprehenr::to_list(for(i in 1:ncol(barcode_counts_tmp)) if(i != ii & i != jj) i)
        
        if (length(new_ix)==0) {
            break
        }

        col_ii <- barcode_counts_tmp[1:nrow(barcode_counts_tmp),ii]
        col_jj <- barcode_counts_tmp[1:nrow(barcode_counts_tmp),jj]  
        
        new_column <- data.matrix(col_ii + col_jj)
        barcode_counts_tmp <- data.matrix(barcode_counts_tmp[1:nrow(barcode_counts_tmp),unlist(new_ix)])
        barcode_counts_new <- cbind(data.matrix(barcode_counts_tmp[1:nrow(barcode_counts_tmp),1:ix-1]),new_column)
        if (ix <= ncol(barcode_counts_tmp)) {
            barcode_counts_new <- cbind(barcode_counts_new, data.matrix(barcode_counts_tmp[1:nrow(barcode_counts_tmp),ix:ncol(barcode_counts_tmp)]))
        }

        node_names_new <- c(node_names[1:ix-1], next_node)
        if (ix <= ncol(barcode_counts_tmp)) {
            node_names_new <- c(node_names_new, node_names[ix:length(node_names)])
        }
           
        next_node <- next_node + 1
        node_names <- node_names_new
        barcode_counts_tmp <- barcode_counts_new
    }
    
    for (i in node_names) {
        parent_map[[i]] <- next_node
    }           
    return(list(parent_map, node_groups, list(X_history, merged_pairs_history, node_names_history)))
}

#' @export
child_map_to_newick_partial <- function(child_map,node) {
    if (length(child_map[[node]]) > 0) {
        recursed <- comprehenr::to_list(for(n in child_map[[node]]) if(TRUE) child_map_to_newick_partial(child_map,n))
        inside <- paste0(recursed, collapse=',')
        return(paste0('(', inside, ')'))
    } else {
        return(node)
    }
}
 
#' @export                           
child_map_to_newick <- function(child_map) {
    root <- toString(length(child_map))
    return(paste0(child_map_to_newick_partial(child_map,root),';'))
}
 
#' @export                                                   
parent_map_to_child_map <- function(parent_map) {
    child_map <- list()
    child_map[[toString(length(parent_map)+1)]] <- list()
    for (n in names(parent_map)) {
        child_map[[n]] <- list()
    }
    for (k in names(parent_map)) {
        child_map[[parent_map[[k]]]] <- c(child_map[[parent_map[[k]]]], k)
    }
    return(child_map)
}                        

#' @export                                                  
plot_hierarchy <- function(parent_map, celltype_names) {
    parent_map_names <- list()
    for (i in 1:length(parent_map)) {
        if ((i <= length(celltype_names)) & (parent_map[[i]] <= length(celltype_names))) {
            parent_map_names[celltype_names[[i]]] <- celltype_names[[parent_map[[i]]]]
        } else if ((i <= length(celltype_names)) & (parent_map[[i]] > length(celltype_names))) {
            parent_map_names[celltype_names[[i]]] <- toString(parent_map[[i]])
        } else {
            parent_map_names[[toString(i)]] <- toString(parent_map[[i]])
        }
    }
    child_map <- parent_map_to_child_map(parent_map_names) 
    newick_tree <- child_map_to_newick(child_map)
    x <- phylogram::read.dendrogram(text = newick_tree)
    pdf(paste0(output_directory,'/neighbor_joining_dendrogram.pdf'))
    plot(x, yaxt = "n")    
    dev.off()
}
                                                      
#' @export                           
mrca <- function(i,j, parent_map) {
    ips <- list(i)
    jps <- list(j)
    while (length(intersect(ips,jps))==0) {
        if (i <= length(parent_map)) {
            ips <- c(ips, list(parent_map[[i]]))
            i <- parent_map[[i]]
        }
        if (j <= length(parent_map)) {
            jps <- c(jps, list(parent_map[[j]]))
            j <- parent_map[[j]]
        }
    }
    return(intersect(ips,jps)[[1]])
}

                                                        
#' @export                            
get_normalized_covariance_boostrap <- function(barcode_counts,N) {
    ccs <- list()
    for (i in 1:N) {
        ix <- ceiling(runif(nrow(barcode_counts),0,nrow(barcode_counts)))
        barcode_counts_tmp <- barcode_counts[ix,1:ncol(barcode_counts)]
        cc <- get_normalized_covariance(barcode_counts_tmp)
        ccs <- c(ccs,list(cc))  
    }
    return(abind::abind(ccs,along=0))
}
                            
#' @export                            
detect_symmetry_violations <- function(barcode_counts, parent_map, symmetry_violation_FDR) {
    N <- 1000
    Xbootstrap <- get_normalized_covariance_boostrap(barcode_counts,N)
    X <- get_normalized_covariance(barcode_counts)
    triples <- list()
    diffs <- list()
    diffs_upper <- list()
    diffs_lower <- list()
    val_pairs <- list()
    for (j in 1:ncol(barcode_counts)) {
        for (i in 1:ncol(barcode_counts)) {
            for (k in 1:ncol(barcode_counts)) {
                if (length(unique(c(i,j,k))) == 3) {
                    p1 <- mrca(i,j, parent_map)
                    p2 <- mrca(i,k, parent_map)
                    if ((p2==mrca(p1,p2,parent_map)) & (p1 != p2)) {
                        triples <- c(triples, list(c(i,j,k)))
                        diffs <- c(diffs, list(X[j,k]-X[i,k]))
                        diffs_upper <- c(diffs_upper, list(quantile(Xbootstrap[1:N,j,k]-Xbootstrap[1:N,i,k], 1-symmetry_violation_FDR)))
                        diffs_lower <- c(diffs_lower, list(quantile(Xbootstrap[1:N,j,k]-Xbootstrap[1:N,i,k], symmetry_violation_FDR)))
                        val_pairs <- c(val_pairs, list(c(X[i,k],X[j,k])))
                    }
                }
            }
        }
    }
    threshold <- median(abs(unlist(diffs)))
    violations <- comprehenr::to_list(for (i in 1:length(diffs_lower)) if(diffs_lower[i]>threshold) triples[[i]])
    print(paste0('Detected ',length(violations),' instances of symmetry violation passing a threshold of ',threshold,' with FDR ', symmetry_violation_FDR))
    return(list(violations,list(unlist(diffs),unlist(diffs_lower),unlist(diffs_upper),threshold,val_pairs)))
}                            

#' @export
plot_violations <- function(output_directory, plot_data) {
    diffs <- plot_data[[1]]
    diffs_lower <- plot_data[[2]]
    diffs_upper <- plot_data[[3]]
    threshold <- plot_data[[4]]
    val_pairs <- plot_data[[5]]
    
    o <- order(diffs)
    o <- unlist(comprehenr::to_list(for(i in o) if(diffs[[i]] > 0) i))
    o1 <- unlist(comprehenr::to_list(for(i in o) if (diffs_lower[[i]] < threshold) i))
    o2 <- unlist(comprehenr::to_list(for(i in o) if (diffs_lower[[i]] >= threshold) i))
        
    pdf(paste0(output_directory,'/symmetry_violations.pdf'))
    par(mfrow=c(1,2), pin=c(2,2))
    plot(1:length(o1), diffs[o1], xlim=c(-1, length(o)), ylim=c(min(diffs_lower[o]),max(diffs_upper[o])), 
         pch=19, col='gray', xlab='Putative symmetric triples', ylab='Diff of normalized covariance')
    par(new=TRUE)
    arrows(x0=1:length(o1), y0=diffs_lower[o1], x1=1:length(o1), y1=diffs_upper[o1], 
           code=3, angle=90, length=0.01, col="gray", lwd=1, xlab='', ylab='')
    par(new=TRUE) 
    arrows(x0=(length(o1)+1):length(o), y0=diffs_lower[o2], x1=(length(o1)+1):length(o), y1=diffs_upper[o2], 
           code=3, angle=90, length=0.01, col="gray", lwd=2, xlab='', ylab='') 
    par(new=TRUE)
    plot((length(o1)+1):length(o), diffs[o2], xlim=c(-1, length(o)), ylim=c(min(diffs_lower[o]),max(diffs_upper[o])), 
         pch=19, col='red', xlab='', ylab='') 
    par(new=TRUE)                    
    plot(c(0,length(o)),c(threshold,threshold), xlim=c(-1, length(o)), ylim=c(min(diffs_lower[o]),max(diffs_upper[o])),
         type="l", col="black", lwd=1, xlab="", ylab="")

    par(new=FALSE, col='gray')
    val_x <- unlist(comprehenr::to_list(for(v in val_pairs) if(TRUE) v[[1]]))
    val_y <- unlist(comprehenr::to_list(for(v in val_pairs) if(TRUE) v[[2]]))
    violations_ix =unlist(comprehenr::to_list(for(i in 1:length(diffs_lower)) if (diffs_lower[[i]] > threshold) i))
    plot(val_x,val_y, pch=19, col='gray', xlim=c(min(val_x), max(val_x)), ylim=c(min(val_y), max(val_y)), 
         xlab='Normalized covariance (j,k)', ylab='Normalized covariance (j,k)')
    par(new=TRUE)                      
    plot(val_x[violations_ix],val_y[violations_ix], xlim=c(min(val_x), max(val_x)), ylim=c(min(val_y), max(val_y)),
         pch=19, col='red', xlab='', ylab='')
    dev.off()

}   

#' @export                                  
tree_dist <- function(i,j,tree_heights, parent_map) {
    k <- mrca(i,j,parent_map)
    if (k %in% c(i,j)) { 
        return(0)
    } else {
        return(tree_heights[[k]])
    }
}                           
                                  
#' @export
get_tree_heights <- function(parent_map, N) {
    tree_heights <- list()
    for (i in 1:N) {
        tree_heights[[i]] <- 0
    }
    for (j in (N+1):(length(parent_map)+1)) {
        prev_heights <- unlist(comprehenr::to_list(for(i in 1:length(parent_map)) if(parent_map[[i]]==j) tree_heights[[i]]))
        tree_heights[[j]] <- max(prev_heights) + 1
    }
    return(tree_heights)
}

#' @export
get_violations <- function(ip,jp,N, parent_map, tree_heights) {
    violations <- list()
    if (mrca(ip,jp,parent_map) %in% c(ip,jp)) {
        return(list())
    }
    for (i in 1:N) {
        for (j in 1:N) {
            for (k in 1:N) {
                if (length(c(i,j,k))==3) {
                    p1 <- mrca(i,j,parent_map)
                    p2 <- mrca(i,k,parent_map)
                    if (p2==mrca(p1,p2,parent_map) & p1 != p2) {
                        if ((0 < tree_dist(k,ip, tree_heights, parent_map)) & (tree_dist(k,ip,tree_heights, parent_map) <  tree_dist(k,jp,tree_heights, parent_map))) {
                            if (tree_dist(i,ip,tree_heights, parent_map)==0 & tree_dist(j,ip,tree_heights, parent_map)>0) {
                                violations <- c(violations, list(c(i,j,k)))
                            }
                        }
                        if (tree_dist(k,ip,tree_heights, parent_map) > tree_dist(k,jp,tree_heights, parent_map)) {
                            if (tree_dist(i,ip,tree_heights, parent_map) > 0 & tree_dist(j,ip,tree_heights, parent_map) == 0) {
                                violations <- c(violations, list(c(i,j,k)))
                            }
                        }
                        if (tree_dist(k,ip,tree_heights, parent_map) == 0) {
                            if (tree_dist(i,ip,tree_heights, parent_map) >  tree_dist(i,jp,tree_heights, parent_map) & tree_dist(i,jp,tree_heights, parent_map) > tree_dist(j,jp,tree_heights, parent_map)) {
                                violations <- c(violations, list(c(i,j,k)))
                            }
                        }
                    }
                }
            }
        }
    }
    return(violations)
}                                       
                                       

#' @export
get_all_node_names <- function(celltype_names, node_groups) {
    all_node_names <- list()
    for (i in 1:length(node_groups)) {
        group_names <- celltype_names[unlist(node_groups[[i]])]
        all_node_names <- c(all_node_names, list(paste0(group_names, collapse=' + ')))
    }
    return(all_node_names)
}   
