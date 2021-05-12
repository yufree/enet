#' Generate metabolites-metabolites network
#' @param meta metabolites dataframe with row as metabolites/peaks and column as samples
#' @param cutoff correlation coefficient cutoff
#' @param m value or vector, specific metabolite or peak to its export network
#' @param name export gml file for Cytoscape
#' @param ... other parameters for `cor`
#' @return list with index of metabolites clusters, metabolites clusters data and edge table
#' @examples
#' data(meta)
#' mm <- getmmnet(meta)
#' @export
getmmnet <- function(meta,cutoff=0.9,m=NULL,name=NULL, ...){
    metacor <- stats::cor(t(meta),...)
    metacor[abs(metacor)<cutoff] <- 0
    df <- data.frame(from=rownames(meta)[which(lower.tri(metacor), arr.ind = T)[, 1]],to=rownames(meta)[which(lower.tri(metacor), arr.ind = T)[, 2]],cor=metacor[lower.tri(metacor)])
    df <- df[abs(df$cor)>0,]
    df$direction <- ifelse(df$cor>0,'positive','negative')
    net <- igraph::graph_from_data_frame(df,directed = F)
    netc <- igraph::components(net)
    message(paste(netc$no, 'metabolites correlation network clusters found'))
    index <- rep(NA,length(rownames(meta)))
    index[match(names(netc$membership),rownames(meta))] <- netc$membership
    message(paste(sum(is.na(index)), 'out of', length(rownames(meta)), 'metabolites have no correlation with other metabolites'))
    re <- list()
    if(is.null(m)){
        for (i in 1:netc$no){
            mcni <- igraph::V(net)$name[netc$membership == i]
            metai <- meta[match(mcni,rownames(meta)),]
            re[[i]] <- metai
        }
    }else{
        n <- unique(netc$membership[igraph::V(net)$name %in% m])
        for(i in 1:length(n)){
            mcni <- igraph::V(net)$name[netc$membership == n[i]]
            metai <- meta[match(mcni,rownames(meta)),]
            re[[i]] <- metai
        }
    }
    if(!is.null(name)){
        igraph::write_graph(net, file=paste0(name,'edgelist.gml'), format = "gml")
    }
    return(list(index=index,data=re,net=df))
}

#' Generate exposures-exposures network
#' @param expo exposures dataframe with row as exposures/peaks and column as samples
#' @param cutoff correlation coefficient cutoff
#' @param nacf NA cutoff for exposure, exposure with NA percentage larger than cutoff will be removed and others will be imputed by knn methods, default 0.2
#' @param e value or vector, specific exposure or peak to its export network
#' @param name export gml file for Cytoscape
#' @param ... other parameters for `cor`
#' @return list with index of exposures clusters, exposures clusters data and edge table
#' @examples
#' data(expo)
#' ee <- geteenet(expo)
#' @export
geteenet <- function(expo,cutoff=0.9,nacf=0.2,e=NULL, name=NULL,...){
    row <- apply(expo,1,function(x)sum(is.na(x))/length(x)<nacf)
    message(paste(sum(row),'exposures could be used for network analysis'))
    expo2 <- expo[row,]
    expo3 <- impute::impute.knn(expo2)$data
    expcor <- stats::cor(t(expo3),...)
    expcor[abs(expcor)<cutoff] <- 0
    df <- data.frame(from=rownames(expo3)[which(lower.tri(expcor), arr.ind = T)[, 1]],to=rownames(expo3)[which(lower.tri(expcor), arr.ind = T)[, 2]],cor=expcor[lower.tri(expcor)])
    df <- df[df$cor>0,]
    if(nrow(df)>0){
        df$direction <- ifelse(df$cor>0,'positive','negative')
        net <- igraph::graph_from_data_frame(df,directed = F)
        netc <- igraph::components(net)
        message(paste(netc$no, 'exposures correlation network clusters found'))
        index <- rep(NA,length(rownames(expo3)))
        index[match(names(netc$membership),rownames(expo3))] <- netc$membership
        message(paste(sum(is.na(index)), 'out of', length(rownames(expo3)), 'exposures have no correlation with other exposures'))
        re <- list()
        if(is.null(e)){
            for (i in 1:netc$no){
                mcni <- igraph::V(net)$name[netc$membership == i]
                expoi <- expo[match(mcni,rownames(expo)),]
                re[[i]] <- expoi
            }
        }else{
            n <- unique(netc$membership[igraph::V(net)$name %in% e])
            for(i in 1:length(n)){
                mcni <- igraph::V(net)$name[netc$membership == n[i]]
                expoi <- expo[match(mcni,rownames(expo)),]
                re[[i]] <- expoi
            }
        }
        if(!is.null(name)){
            igraph::write_graph(net, file=paste0(name,'edgelist.gml'), format = "gml")
        }
        return(list(index=index,data=re,net=df))
    }else{
        message('exposures have no correlation with each others.')
        return(NULL)
        }
}

#' Generate metabolites-exposures network
#' @param meta metabolites dataframe with row as metabolites/peaks and column as samples
#' @param expo exposures dataframe with row as exposures/peaks and column as samples
#' @param cutoffm correlation coefficient cutoff for metabolites-metabolites network
#' @param cutoffe correlation coefficient cutoff for exposures-exposures network
#' @param nacf NA cutoff for exposure, exposure with NA percentage larger than cutoff will be removed and others will be imputed by knn methods, default 0.2
#' @param ... other parameters for `cor`
#' @return list with table of metabolites-exposures, metabolites-exposures network clusters data and edge table
#' @examples
#' data(expo)
#' data(meta)
#' me <- getemnet(meta,expo)
#' @export
getemnet <- function(meta,expo,cutoffm=0.9,cutoffe=0.6,nacf=0.2,...){
    colindex <- match(colnames(expo),colnames(meta))
    colindex2 <- match(colnames(meta),colnames(expo))
    meta <- meta[,colindex[!is.na(colindex)]]
    expo <- expo[,colindex2[!is.na(colindex2)]]
    tst <- function(n, ...) ...elt(n)
    m <- tst(1,getmmnet(meta,cutoff = cutoffm,...),geteenet(expo,cutoff = cutoffe,nacf=nacf,...))
    e <- tst(2,getmmnet(meta,cutoff = cutoffm,...),geteenet(expo,cutoff = cutoffe,nacf=nacf,...))
    df <- m$net
    dfe <- e$net

    ex <- ic <- dt <- c()
    for(i in 1:nrow(expo)){
        fit1 <- limma::lmFit(log2(meta+1), stats::model.matrix(~expo[i,]))
        fit2  <- limma::eBayes(fit1)

        idx <- limma::topTable(fit2,coef=2,adjust.method = "BH",sort.by = 'none',number = nrow(fit2))

        ic <- c(ic,rownames(meta)[idx$adj.P.Val<0.05])
        ex <- c(ex,rep(rownames(expo)[i],sum(idx$adj.P.Val<0.05)))
        dt <- c(dt,ifelse(idx$t[idx$adj.P.Val<0.05]>0,'positive','negative'))
    }
    dfme <- cbind.data.frame(from=ic,to=ex,cor=1,direction=dt)
    # for individual exposure
    tab <- table(dfme$from,dfme$to)
    li <- apply(tab,1,sum)
    message(paste(length(li[li>1]),'peaks with multiple exposures associated'))
    for(i in unique(li)){
        message(paste(length(li[li==i]),'peaks with', i ,'exposures associated'))
    }
    # li <- li[li>2]
    net <- igraph::graph_from_data_frame(df,directed = F)
    netc <- igraph::components(net)
    lire <- list()
    for (i in 1:length(li)){
        dfmei <- dfme[dfme$from %in% names(li[i]),]
        if(sum(names(netc$membership) %in% names(li[i]))!=0){
            c <- names(netc$membership[netc$membership == netc$membership[names(netc$membership) %in% names(li[i])]])
            dfi <- df[df$from %in% c|df$to %in% c,]
            re <- rbind.data.frame(dfi,dfmei)
        }else{
            dfmei$m <- 'none'
            re <- dfmei
        }
        lire[[i]] <- re
    }
    message(paste(sum(sapply(lire, dim)[2,]==4),'peaks connected with other peaks and exposures'))
    message(paste(sum(sapply(lire, dim)[2,]==5),'peaks connected with exposures only'))

    # for exposure cluster
    lire2 <- list()
    if(!is.null(dfe)){
        ec <- rownames(expo)[!is.na(e$index)]
        ecc <- e$index[!is.na(e$index)]
        mc <- rownames(meta)[!is.na(m$index)]
        mcc <- m$index[!is.na(m$index)]
        if(length(unique(ecc))>1){
            for(i in 1:length(unique(ecc))){
                eci <- ec[ecc==i]
                dfmei <- dfme[dfme$to %in% eci,]
                mcci <- mcc[mc %in% dfmei$from]
                dfmmi <- df[df$from %in% mc[mcc %in% unique(mcci)]|df$to %in% dfmei$mc[mcc %in% unique(mcci)],]
                dfmeci <- rbind.data.frame(dfmmi,dfmei)
                lire2[[i]] <- dfmeci
            }
        }else{
            dfmei <- dfme[dfme$to %in% ec,]
            mcci <- mcc[mc %in% dfmei$from]
            dfmmi <- df[df$from %in% mc[mcc %in% unique(mcci)]|df$to %in% dfmei$mc[mcc %in% unique(mcci)],]
            dfmeci <- rbind.data.frame(dfmmi,dfmei)
            lire2 <- dfmeci
        }
        li <- list(me=tab,data=lire,cluster = lire2, metaexp=dfme)
    }else{
        li <- list(me=tab,data=lire, metaexp=dfme)
    }
    return(li)
}

#' Gatekeeper discovery
#' @param meta metabolites dataframe with row as metabolites/peaks and column as samples
#' @param expo exposures dataframe with row as exposures/peaks and column as samples
#' @param cutoff correlation coefficient cutoff for metabolites-metabolites network
#' @param multiple logical, only output gatekeeper with multiple exposures
#' @param ... other parameters for `cor`
#' @return list with table of metabolites-exposures, metabolites-exposures network clusters data and edge table
#' @examples
#' data(expo)
#' data(meta)
#' gk <- getgk(meta,expo)
#' @export
getgk <- function(meta,expo,cutoff=0.9,multiple=FALSE,...){
    colindex <- match(colnames(expo),colnames(meta))
    colindex2 <- match(colnames(meta),colnames(expo))
    meta <- meta[,colindex[!is.na(colindex)]]
    expo <- expo[,colindex2[!is.na(colindex2)]]

    m <- getmmnet(meta,cutoff = cutoff,...)
    df <- m$net

    ex <- ic <- dt <- c()
    for(i in 1:nrow(expo)){
        fit1 <- limma::lmFit(log2(meta+1), stats::model.matrix(~expo[i,]))
        fit2  <- limma::eBayes(fit1)

        idx <- limma::topTable(fit2,coef=2,adjust.method = "BH",sort.by = 'none',number = nrow(fit2))

        ic <- c(ic,rownames(meta)[idx$adj.P.Val<0.05])
        ex <- c(ex,rep(rownames(expo)[i],sum(idx$adj.P.Val<0.05)))
        dt <- c(dt,ifelse(idx$t[idx$adj.P.Val<0.05]>0,'positive','negative'))
    }
    dfme <- cbind.data.frame(from=ic,to=ex,cor=1,direction=dt)
    # for individual exposure
    tab <- table(dfme$from,dfme$to)
    li <- apply(tab,1,sum)
    if(multiple){
        li <- li[li>1]
    }

    net <- igraph::graph_from_data_frame(df,directed = F)
    netc <- igraph::components(net)
    lire <- list()
    for (i in 1:length(li)){
        dfmei <- dfme[dfme$from %in% names(li[i]),]
        if(sum(names(netc$membership) %in% names(li[i]))!=0){
            c <- names(netc$membership[netc$membership == netc$membership[names(netc$membership) %in% names(li[i])]])
            dfi <- df[df$from %in% c|df$to %in% c,]
            re <- rbind.data.frame(dfi,dfmei)
        }else{
            dfmei$m <- 'none'
            re <- dfmei
        }
        lire[[i]] <- re
    }
    lire2 <- lire[sapply(lire, dim)[2,]==4]
    liname <- names(li)[sapply(lire, dim)[2,]==4]
    dfx <- dfme[dfme$from %in% liname,]
    tab <- table(dfx$from,dfx$to)
    dfme2 <- do.call(rbind.data.frame,lire2)
    dfme2 <- dfme2[!duplicated(dfme2),]
    name <- unique(c(dfme2$from,dfme2$to)[!c(dfme2$from,dfme2$to)%in%dfme$to])
    message(paste(length(name), 'peaks were involved.'))
    message(paste(sum(sapply(lire, dim)[2,]==4),'peaks could be gatekeepers.'))
    li <- list(me=tab,data=lire2, metaexp=dfme2)
    return(li)
}