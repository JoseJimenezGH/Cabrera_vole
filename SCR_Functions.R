
SCRgof<-function (out, nx = 6, ny = 6, traplocs = NULL, buffer = 2, Xl = NULL,
    Xu = NULL, Yl = NULL, Yu = NULL)
{
    S <- out$s
    Sxout <- S[, , 1]
    Syout <- S[, , 2]
    z <- out$w
    Xl <- min(traplocs[, 1]) - buffer
    Xu <- max(traplocs[, 1]) + buffer
    Yl <- min(traplocs[, 2]) - buffer
    Yu <- max(traplocs[, 2]) + buffer
    niter <- nrow(z)
    xg <- seq(Xl, Xu, , nx)
    yg <- seq(Yl, Yu, , ny)
    area <- (Xu - Xl) * (Yu - Yl)
    Sxout2 <- cut(Sxout[z == 1], breaks = xg)
    Syout2 <- cut(Syout[z == 1], breaks = yg)
    Dn <- table(Sxout2, Syout2)/niter
    op <- par(mar = c(4, 4, 4, 6))
    on.exit(par(op))
    image(xg, yg, Dn, col = terrain.colors(10))
    image.scale(Dn, col = terrain.colors(10))
    title("Tamano de la poblacion local \n(individuos por unidad de superficie)")
    stat2 <- statsim2 <- stat <- statsim <- rep(NA, niter)
    for (i in 1:niter) {
        inside <- (Sxout[i, ] < Xu) & (Sxout[i, ] > Xl) & (Syout[i,
            ] < Yu) & (Syout[i, ] > Yl)
        Dn <- table(cut(Sxout[i, ][z[i, ] == 1 & inside], breaks = xg),
            cut(Syout[i, ][z[i, ] == 1 & inside], breaks = yg))
        Dnv <- Dn[1:length(Dn)]
        E <- mean(Dnv)
        stat[i] <- (var(Dnv)/mean(Dnv))
        stat2[i] <- sum((sqrt(Dnv) - sqrt(E))^2)
        Sxsim <- runif(sum(z[i, ][inside]), Xl, Xu)
        Sysim <- runif(sum(z[i, ][inside]), Yl, Yu)
        Dnsim <- table(cut(Sxsim, breaks = xg), cut(Sysim, breaks = yg))
        Dnsimv <- Dnsim[1:length(Dnsim)]
        statsim[i] <- (var(Dnsimv)/mean(Dnsimv))
        statsim2[i] <- sum((sqrt(Dnsimv) - sqrt(E))^2)
    }
    out <- cbind(data = stat2, newdata = statsim2)
    cat("Indice de dispersion observado: ", mean(stat), fill = TRUE)
    cat("Indice de dispersion simulado: ", mean(statsim), fill = TRUE)
    cat("P-valor del Indice of dispersion: ", mean(statsim > stat),
        fill = TRUE)
    cat("P-valor(2) freeman-tukey: ", mean(statsim2 > stat2), fill = TRUE)
    invisible(out)
}

pGauss1<-function (parms, Dmat) 
{
    a0 <- parms[1]
    sigma <- parms[2]
    p <- plogis(parms[1]) * exp(-(1/(2 * parms[2] * parms[2])) * 
        Dmat * Dmat)
    p
}


e2dist<-function (x, y)
{
    i <- sort(rep(1:nrow(y), nrow(x)))
    dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
    matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}


image.scale<-function (z, col, x, y = NULL, size = NULL, digits = 2, cex=0.8,
  labels = c("breaks", "ranges"))
{
    n <- length(col)
    usr <- par("usr")
    mx <- mean(usr[1:2])
    my <- mean(usr[3:4])
    dx <- diff(usr[1:2])
    dy <- diff(usr[3:4])
    if (missing(x))
        x <- mx + 1.05 * dx/2
    else if (is.list(x)) {
        if (length(x$x) == 2)
            size <- c(diff(x$x), -diff(x$y)/n)
        y <- x$y[1]
        x <- x$x[1]
    }
    else x <- x[1]
    if (is.null(size))
        if (is.null(y)) {
            size <- 0.618 * dy/n
            y <- my + 0.618 * dy/2
        }
        else size <- (y - my) * 2/n
    if (length(size) == 1)
        size <- rep(size, 2)
    if (is.null(y))
        y <- my + n * size[2]/2
    i <- seq(along = col)
    rect(x, y - i * size[2], x + size[1], y - (i - 1) * size[2],
        col = rev(col), xpd = TRUE)
    rng <- range(z, na.rm = TRUE)
    bks <- seq(from = rng[2], to = rng[1], length = n + 1)
    bks <- formatC(bks, format = "f", digits = digits)
    labels <- match.arg(labels)
    if (labels == "breaks")
        ypts <- y - c(0, i) * size[2]
    else {
        bks <- paste(bks[-1], bks[-(n + 1)], sep = " - ")
        ypts <- y - (i - 0.5) * size[2]
    }
    text(x = x + 1.2 * size[1], y = ypts, labels = bks, 
      adj = ifelse(size[1] > 0, 0, 1), xpd = TRUE, cex=0.8)
}


SCRdensity<-function (obj, nx = 30, ny = 30, Xl = NULL, Xu = NULL, Yl = NULL,
    Yu = NULL, scalein = 100, scaleout = 100, ncolors = 10)
{
    Sxout <- obj$Sx
    Syout <- obj$Sy
    z <- obj$z
    niter <- nrow(z)
    if (is.null(Xl)) {
        Xl <- min(Sxout) * 0.999
        Xu <- max(Sxout) * 1.001
        Yl <- min(Syout) * 0.999
        Yu <- max(Syout) * 1.001
    }
    xg <- seq(Xl, Xu, , nx)
    yg <- seq(Yl, Yu, , ny)
    Sxout <- cut(Sxout[z == 1], breaks = xg)
    Syout <- cut(Syout[z == 1], breaks = yg)
    Dn <- table(Sxout, Syout)/niter
    area <- (yg[2] - yg[1]) * (xg[2] - xg[1]) * scalein
    Dn <- (Dn/area) * scaleout
    cat("mean: ", mean(Dn), fill = TRUE)
    par(mar = c(3, 3, 3, 6))
    image(xg, yg, Dn, col = terrain.colors(ncolors))
    image.scale(Dn, col = terrain.colors(ncolors),cex=0.75)
    box()
    return(list(grid = cbind(xg, yg), Dn = Dn))
}


SCRdensityJJ1<-function (obj, nx = 30, ny = 30, Xl = NULL, Xu = NULL, Yl = NULL,
    Yu = NULL, scalein = 100, scaleout = 100,  ncolors = 10,
    col=terrain.colors(ncolors), asp=1)
{
    Sxout <- obj$Sx
    Syout <- obj$Sy
    z <- obj$z
    niter <- nrow(z)
    if (is.null(Xl)) {
        Xl <- min(Sxout) * 0.999
        Xu <- max(Sxout) * 1.001
        Yl <- min(Syout) * 0.999
        Yu <- max(Syout) * 1.001
    }
    xg <- seq(Xl, Xu, , nx)
    yg <- seq(Yl, Yu, , ny)
    Sxout <- cut(Sxout[z == 1], breaks = xg)
    Syout <- cut(Syout[z == 1], breaks = yg)
    Dn <- table(Sxout, Syout)/niter
    area <- (yg[2] - yg[1]) * (xg[2] - xg[1]) * scalein
    Dn <- (Dn/area) * scaleout
    cat("mean: ", mean(Dn), fill = TRUE)
    par(mar = c(3, 3, 3, 6))
    image(xg, yg, Dn, col=terrain.colors(ncolors),asp=1)
    image.scale(Dn, col=terrain.colors(ncolors))
    box()
    return(list(xg=xg, yg=yg, grid = cbind(xg, yg), Dn = Dn))
}


SCRdensityJJ2 <- function (obj, nx = 30, ny = 30, Xl = NULL, Xu = NULL, Yl = NULL,
    Yu = NULL, scalein = 100, scaleout = 100, ncolors = 10, asp=1)
{
    Sxout <- obj$Sx
    Syout <- obj$Sy
    z <- obj$z
    niter <- nrow(z)
    if (is.null(Xl)) {
        Xl <- min(Sxout) * 0.999
        Xu <- max(Sxout) * 1.001
        Yl <- min(Syout) * 0.999
        Yu <- max(Syout) * 1.001
    }
    xg <- seq(Xl, Xu, , nx)
    yg <- seq(Yl, Yu, , ny)
    Sxout <- cut(Sxout[z == 1], breaks = xg)
    Syout <- cut(Syout[z == 1], breaks = yg)
    Dn <- table(Sxout, Syout)/niter
    area <- (yg[2] - yg[1]) * (xg[2] - xg[1]) * scalein
    Dn <- (Dn/area) * scaleout
    cat("mean: ", mean(Dn), fill = TRUE)
    par(mar = c(3, 3, 3, 6))
    image(xg, yg, Dn, col=gray.colors(ncolors, start=1, end=0),asp=1)
    image.scale(Dn, col=gray.colors(ncolors, start=1, end=0))
    box()
    return(list(grid = cbind(xg, yg), Dn = Dn))
}



hra<-function (func, parms, plot = TRUE, xlim, ylim, ng = 100, target.area = NULL, 
    tol = 0.001) 
{
    s <- c((xlim[2] - xlim[1])/2, (ylim[2] - ylim[1])/2)
    x1 <- rep(seq(xlim[1], xlim[2], , ng), ng)
    x2 <- sort(rep(seq(ylim[1], ylim[2], , ng), ng))
    delta <- min(diff(x1[1:10]))
    x1 <- rep(seq(xlim[1] - delta/2, xlim[2] + delta/2, delta), 
        ng)
    x2 <- sort(rep(seq(ylim[1] - delta/2, ylim[2] + delta/2, 
        delta), ng))
    X <- cbind(x1, x2)
    D <- sqrt((s[1] - x1)^2 + (s[2] - x2)^2)
    p <- func(parms, D)
    if (plot) {
        spatial.plot(X, p)
    }
    psi <- p/sum(p)
    if (is.null(target.area)) {
        x0 <- 0.2
        repeat {
            in.hr <- D <= x0
            total <- sum(psi[in.hr])
            if (total >= 0.95) {
                print(x0)
                break
            }
            x0 <- x0 * (1 + tol)
        }
        radius <- x0
        cat("radio para alcanzar el 95% del area: ", radius, fill = TRUE)
        area <- pi * radius^2
        cat("Area de campeo: ", area, fill = TRUE)
        return(area)
    }
    if (!is.null(target.area)) {
        if (is.null(target.area)) {
            cat("need target.area", fill = TRUE)
            goose.egg <- NULL
            return(goose.egg)
        }
        obj <- function(beta2) {
            p <- func(c(parms[1], beta2), D)
            psi <- p/sum(p)
            x0 <- 0.1
            repeat {
                in.hr <- D <= x0
                total <- sum(psi[in.hr])
                if (total >= 0.95) {
                  break
                }
                x0 <- x0 * (1 + tol)
            }
            hr.area <- pi * x0 * x0
            ss <- (hr.area - target.area)^2
            ss
        }
        tmp <- optimize(obj, interval = c(0.01, 5))
        beta2 <- tmp$minimum
        cat("Valor del parm[2] para alcanzar el 95% del area de campeo ", 
            target.area, " : ", beta2, fill = TRUE)
        return(beta2)
    }
}

spatial.plot<- function (x, y, add = FALSE, cx = 1, col = "gray") 
{
    nc <- as.numeric(cut(y, 10))
    if (!add) 
        plot(x, pch = " ", asp = 1)
    if (col == "gray") {
        cc <- seq(3, 17, , 10)/20
        cc <- gray(cc)
    }
    else cc <- terrain.colors(10)
    points(x, pch = 20, col = cc[nc], cex = cx)
    image.scale(y, col = cc)
}


SCR23darray<-function (edf, tdf)
{
    nind <- max(edf[, 2])
    ntraps <- nrow(tdf)
    nperiods <- ncol(tdf) - 3
    per.id <- as.numeric(dimnames(tdf)[[2]][4:ncol(tdf)])
    ind.id <- edf[, 2]
    trap.id <- edf[, 4]
    if (length(per.id) != length(min(per.id):max(per.id))) {
        x <- 1:nperiods
        names(x) <- as.character(per.id)
        per.id <- x[as.character(edf[, 3])]
    }
    else {
        per.id <- edf[, 3]
    }
    y <- array(0, c(nind, ntraps, nperiods))
    tmp <- cbind(ind.id, trap.id, per.id)
    y[tmp] <- 1
    y
}


SCRsmy<-function (y3d)
{
    nind <- dim(y3d)[1]
    totcaps <- nperiods <- sprecaps <- rep(NA, nind)
    for (i in 1:nind) {
        x <- y3d[i, , ]
        ntraps <- sum(apply(x, 2, sum) > 0)
        ncaps <- sum(x)
        nperiods[i] <- sum(apply(x, 1, sum) > 0)
        sprecaps[i] <- ifelse(ntraps > 1, 1, 0) * ncaps
        totcaps[i] <- sum(x)
    }
    cat("Total de capturas: ", sum(totcaps), fill = TRUE)
    cat("Recapturas espaciales: ", sum(sprecaps), fill = TRUE)
    cat("Eventos de captura ordinarios: ", sum(nperiods), fill = TRUE)
    cat("Capturas perdidas en el modelo no espacial: ", sum(totcaps) -
        sum(nperiods), fill = TRUE)
}


make.grid<-function (ll = NA, minx = NA, maxx = NA, miny = NA, maxy = NA,
nx = 40, ny = NULL, buffer = 0)
  {
    if (is.null(ny))
    ny <- nx
    if (!is.na(ll)) {
    minx <- min(ll[, 1])
    maxx <- max(ll[, 1])
    miny <- min(ll[, 2])
    maxy <- max(ll[, 2])
    bx <- (maxx - minx) * buffer
    by <- (maxy - miny) * buffer
    minx <- minx - bx
    maxx <- maxx + bx
    miny <- miny - by
    maxy <- maxy + by
  }
  x <- sort(rep(seq(minx, maxx, , nx), ny))
  y <- rep(seq(maxy, miny, , ny), nx)
  cbind(x, y)
}


rot<-function (m) 
{
    nr <- nrow(m)
    nc <- ncol(m)
    v <- matrix(NA, nrow = nc, ncol = nr)
    for (i in 1:nr) {
        v[, nr - (i - 1)] <- m[i, ]
    }
    v
}

make.scrFrame<-function (caphist, traps, indCovs = NULL, trapCovs = NULL, sigCovs = NULL,
    trapOperation = NULL, telemetry = NULL, rsfDF = NULL, type = "scr")
{
    if (any(is.null(caphist), is.null(traps)))
        stop("caphist and trap must be provided")
    if (!is.list(caphist))
        stop("caphist must be a list")
    n.sessions <- length(caphist)
    caphist.dimensions <- sapply(caphist, dim)
    if (nrow(caphist.dimensions) == 2)
        caphist.dimensions <- rbind(caphist.dimensions, 1)
    for (i in 1:n.sessions) {
        caphist[[i]] <- array(caphist[[i]], dim = caphist.dimensions[,
            i])
        all.zero <- apply(apply(caphist[[i]], c(1, 3), sum),
            1, sum)
        if (any(all.zero == 0)) {
            cat("At least one individual has an all-zero encounter history",
                fill = TRUE)
            cat("Make sure this is ok...", fill = TRUE)
        }
    }
    if (!is.null(indCovs)) {
        if (!is.list(indCovs))
            stop("indCovs must be a list")
        if (any(!sapply(indCovs, is.data.frame)))
            stop("indCovs must be a list of dataframes")
        if (length(indCovs) != length(caphist))
            stop("number of sessions in indCovs does not match caphist")
        check.dim <- sapply(indCovs, nrow)
        if (any(check.dim != caphist.dimensions[1, ]))
            stop("number of individuals in indCovs does not match caphist")
        if (!("rmv" %in% indCovs[[1]])) {
            for (i in 1:length(indCovs)) {
                indCovs[[i]]$removed <- dim(caphist[[i]])[3]
            }
        }
    }
    else {
        indCovs <- list()
        for (i in 1:length(caphist)) {
            indCovs[[i]] <- data.frame(removed = rep(dim(caphist[[i]])[3],
                dim(caphist[[i]])[1]))
        }
    }
    if (!is.list(traps))
        stop("traps must be a list")
    if (length(traps) != length(caphist))
        stop("number of sessions in traps does not match caphist")
    check.dim <- sapply(traps, nrow)
    if (!all(check.dim == caphist.dimensions[2, ]))
        stop("number of traps does not match caphist")
    if (!is.null(trapCovs)) {
        if (!is.list(trapCovs))
            stop("trapCovs must be a list")
        if (any(!sapply(trapCovs, is.list)))
            stop("trapCovs must be a list of lists")
        if (any(!unlist(sapply(trapCovs, function(x) sapply(x,
            is.data.frame)))))
            stop("trapCovs must be a list of dataframes")
        if (length(trapCovs) != length(caphist))
            stop("number of sessions in trapCovs does not match caphist")
        check.dim <- lapply(trapCovs, function(x) sapply(x, nrow))
        for (i in 1:length(check.dim)) {
            if (!all(check.dim[[i]] == caphist.dimensions[2,
                i]))
                stop("number of traps does not match caphist")
        }
    }
    if (!is.null(sigCovs)) {
        if (nrow(sigCovs) != length(caphist))
            stop("number of rows in sigCovs does not match number of sessions")
        if (!"session" %in% colnames(sigCovs)) {
            sigCovs$session <- factor(1:n.sessions)
        }
        if (!is.null(indCovs)) {
            if ("sex" %in% colnames(indCovs[[1]])) {
                sigCovs <- sigCovs[rep(1:n.sessions, 2), , drop = F]
                rownames(sigCovs) <- NULL
                sigCovs$sex <- factor(rep(c("female", "male"),
                  each = n.sessions))
            }
        }
    }
    else {
        sigCovs <- data.frame(session = factor(1:n.sessions))
        if (!is.null(indCovs)) {
            if ("sex" %in% colnames(indCovs[[1]])) {
                sigCovs <- sigCovs[rep(1:n.sessions, 2), , drop = F]
                rownames(sigCovs) <- NULL
                sigCovs <- sigCovs[rep(1:n.sessions, 2), , drop = F]
                rownames(sigCovs) <- NULL
                sigCovs$sex <- factor(rep(c("female", "male"),
                  each = n.sessions))
            }
        }
    }
    if (!is.null(trapOperation)) {
        if (!is.list(trapOperation))
            stop("trapOperation must be a list")
        if (length(trapOperation) != length(caphist))
            stop("number of sessions in trapOperation does not match caphist")
        check.dim <- sapply(trapOperation, nrow)
        if (!all(check.dim == caphist.dimensions[2, ]))
            stop("number of traps does not match caphist")
    }
    max.dist <- NULL
    for (i in 1:length(caphist)) {
        for (j in 1:nrow(caphist[[i]])) {
            if (dim(caphist[[i]])[3] > 1) {
                where <- apply(caphist[[i]][j, , ], 1, sum) >
                  0
            }
            else {
                where <- caphist[[i]][j, , ] > 0
            }
            if (sum(where) > 1)
                max.dist <- c(max.dist, max(0, dist(traps[[i]][where,
                  c("X", "Y")]), na.rm = T))
        }
    }
    mmdm <- mean(max.dist[max.dist > 0], na.rm = T)
    mdm <- max(max.dist, na.rm = T)
    if (!is.null(telemetry)) {
        if (!is.list(telemetry$fixfreq))
            stop("telemetry$fixfreq must be a list")
        fixfreq.dimensions <- sapply(telemetry$fixfreq, dim)
        if (nrow(fixfreq.dimensions) == 2)
            fixfreq.dimensions <- rbind(fixfreq.dimensions, 1)
        if (!is.null(telemetry$indCovs)) {
            if (!is.list(telemetry$indCovs))
                stop("telemetry$indCovs must be a list")
            if (any(!sapply(telemetry$indCovs, is.data.frame)))
                stop("telemetry$indCovs must be a list of dataframes")
            if (length(telemetry$indCovs) != length(telemetry$fixfreq))
                stop("number of sessions in telemetry$indCovs does not match telemetry$fixfreq")
            check.dim <- sapply(telemetry$indCovs, nrow)
            if (any(check.dim != fixfreq.dimensions[1, ]))
                stop("number of individuals in telemetry$indCovs does not match telemetry$fixfreq")
            if (any(!names(indCovs[[1]]) %in% c(names(telemetry$indCovs[[1]]),
                "removed")))
                stop("indCovs do not match between capture and telemetry data")
        }
        if (!is.null(telemetry$cap.tel)) {
            if (!is.list(telemetry$cap.tel))
                stop("telemetry$indCovs must be a list")
            warning("make sure captured individuals w/ collars sorted first!")
        }
    }
    if (!is.null(rsfDF)) {
        library(FNN)
        rsfCovs <- names(rsfDF[[1]][, -c(1:2), drop = F])
        if (is.null(trapCovs)) {
            trapCovs <- list()
            length(trapCovs) <- n.sessions
            for (s in 1:n.sessions) {
                trap.grid <- as.vector(get.knnx(rsfDF[[s]][,
                  c("X", "Y")], traps[[s]][, c("X",
                  "Y")], 1)$nn.index)
                trapCovs[[s]] <- list()
                length(trapCovs[[s]]) <- caphist.dimensions[3,
                  s]
                for (k in 1:caphist.dimensions[3, s]) {
                  trapCovs[[s]][[k]] <- data.frame(rsfDF[[s]][trap.grid,
                    rsfCovs])
                  names(trapCovs[[s]][[k]]) <- rsfCovs
                }
            }
        }
        else {
            for (s in 1:n.sessions) {
                if (any(!rsfCovs %in% trapCovs[[s]][[1]])) {
                  miss.rsfCovs <- rsfCovs[which(!rsfCovs %in%
                    trapCovs[[s]][[1]])]
                  trap.grid <- as.vector(get.knnx(rsfDF[[s]][,
                    c("X", "Y")], traps[[s]][, c("X",
                    "Y")], 1)$nn.index)
                  for (k in 1:caphist.dimensions[3, s]) {
                    newtrapCovs <- data.frame(rsfDF[[s]][trap.grid,
                      miss.rsfCovs])
                    names(newtrapCovs) <- miss.rsfCovs
                    trapCovs[[s]][[k]] <- data.frame(trapCovs[[s]][[k]],
                      newtrapCovs)
                  }
                }
            }
        }
    }
    scrFrame <- list(caphist = caphist, traps = traps, indCovs = indCovs,
        trapCovs = trapCovs, sigCovs = sigCovs, trapOperation = trapOperation,
        occasions = caphist.dimensions[3, ], type = type, mmdm = mmdm,
        mdm = mdm, telemetry = telemetry)
    class(scrFrame) <- "scrFrame"
    return(scrFrame)
}



# Modified from scrbook and oSCR packages
spiderplotJJ<-function (y, X, ax = TRUE, buffer=buffer, lwd=1)
{
   traps.df<-data.frame(X=X[,1],Y=X[,2])
   scrFrame  <- make.scrFrame(caphist=list(y),
                             traps=list(traps.df),
                             trapCovs=NULL,
                             trapOperation=NULL)
    class(scrFrame) <- "scrFrame"
    op <- par(no.readonly = TRUE)
    all.ind.xy <- list()
    mean.loc <- list()
    mu.x <- list()
    mu.y <- list()
    for (s in 1:length(scrFrame$caphist)) {
        s.ind.xy <- NULL
        s.ind <- NULL
        tmp.ch <- scrFrame$caphist[[s]]
        tmp.tr <- scrFrame$traps[[s]]
        for (i in 1:nrow(tmp.ch)) {
            if (dim(tmp.ch)[3] > 1) {
                pick <- apply(tmp.ch[i, , ], 1, sum) > 0
            }
            else {
                pick <- tmp.ch[i, , ] > 0
            }
            s.ind.xy <- rbind(s.ind.xy, tmp.tr[rep(1:nrow(tmp.tr),
                pick), c("X", "Y")])
            s.ind <- c(s.ind, rep(i, sum(pick)))
        }
        all.ind.xy[[s]] <- data.frame(ind = s.ind, x = s.ind.xy[,
            1], y = s.ind.xy[, 2])
        mu.x[[s]] <- tapply(all.ind.xy[[s]]$x, all.ind.xy[[s]]$ind,
            mean)
        mu.y[[s]] <- tapply(all.ind.xy[[s]]$y, all.ind.xy[[s]]$ind,
            mean)
    }
    par(oma = c(0, 0, 0, 0))
    for (s in 1:length(scrFrame$caphist)) {
        Xl<-min(X[,1])-buffer
        Xu<-max(X[,1])+buffer
        Yl<-min(X[,2])-buffer
        Yu<-max(X[,2])+buffer
        xlim<-c(Xl,Xu)
        ylim<-c(Yl,Yu)

        #plot(scrFrame$traps[[s]][, c("X", "Y")], asp = 1, type = "n",
        #    las = 1, axes = ax, xlab = "", ylab = "", xlim=xlim, ylim=ylim)
        clr <- sample(colors(), nrow(tmp.ch+20))
        box(bty = "o")
        for (i in 1:nrow(tmp.ch)) {
            to.x <- all.ind.xy[[s]]$x[all.ind.xy[[s]]$ind %in%
                i]
            to.y <- all.ind.xy[[s]]$y[all.ind.xy[[s]]$ind %in%
                i]
            segments(rep(mu.x[[s]][i], length(to.x)), rep(mu.y[[s]][i],
                length(to.y)), to.x, to.y, col = "black", lwd = lwd)
        }
        #points(scrFrame$traps[[s]][, c("X", "Y")], pch = "+",
        #    cex = 1)                            # Tamaño  # fondo  # grosor
        points(mu.x[[s]], mu.y[[s]], pch = 21, cex = 1.75, bg = clr, lwd=2.25)
    }
    par(op)
}


make.grid<-function (ll = NA, minx = NA, maxx = NA, miny = NA, maxy = NA,
nx = 40, ny = NULL, buffer = 0)
  {
    if (is.null(ny))
    ny <- nx
    if (!is.na(ll)) {
    minx <- min(ll[, 1])
    maxx <- max(ll[, 1])
    miny <- min(ll[, 2])
    maxy <- max(ll[, 2])
    bx <- (maxx - minx) * buffer
    by <- (maxy - miny) * buffer
    minx <- minx - bx
    maxx <- maxx + bx
    miny <- miny - by
    maxy <- maxy + by
  }
  x <- sort(rep(seq(minx, maxx, , nx), ny))
  y <- rep(seq(maxy, miny, , ny), nx)
  cbind(x, y)
}


rot<-function (m) 
{
    nr <- nrow(m)
    nc <- ncol(m)
    v <- matrix(NA, nrow = nc, ncol = nr)
    for (i in 1:nr) {
        v[, nr - (i - 1)] <- m[i, ]
    }
    v
}

make.scrFrame<-function (caphist, traps, indCovs = NULL, trapCovs = NULL, sigCovs = NULL,
    trapOperation = NULL, telemetry = NULL, rsfDF = NULL, type = "scr")
{
    if (any(is.null(caphist), is.null(traps)))
        stop("caphist and trap must be provided")
    if (!is.list(caphist))
        stop("caphist must be a list")
    n.sessions <- length(caphist)
    caphist.dimensions <- sapply(caphist, dim)
    if (nrow(caphist.dimensions) == 2)
        caphist.dimensions <- rbind(caphist.dimensions, 1)
    for (i in 1:n.sessions) {
        caphist[[i]] <- array(caphist[[i]], dim = caphist.dimensions[,
            i])
        all.zero <- apply(apply(caphist[[i]], c(1, 3), sum),
            1, sum)
        if (any(all.zero == 0)) {
            cat("At least one individual has an all-zero encounter history",
                fill = TRUE)
            cat("Make sure this is ok...", fill = TRUE)
        }
    }
    if (!is.null(indCovs)) {
        if (!is.list(indCovs))
            stop("indCovs must be a list")
        if (any(!sapply(indCovs, is.data.frame)))
            stop("indCovs must be a list of dataframes")
        if (length(indCovs) != length(caphist))
            stop("number of sessions in indCovs does not match caphist")
        check.dim <- sapply(indCovs, nrow)
        if (any(check.dim != caphist.dimensions[1, ]))
            stop("number of individuals in indCovs does not match caphist")
        if (!("rmv" %in% indCovs[[1]])) {
            for (i in 1:length(indCovs)) {
                indCovs[[i]]$removed <- dim(caphist[[i]])[3]
            }
        }
    }
    else {
        indCovs <- list()
        for (i in 1:length(caphist)) {
            indCovs[[i]] <- data.frame(removed = rep(dim(caphist[[i]])[3],
                dim(caphist[[i]])[1]))
        }
    }
    if (!is.list(traps))
        stop("traps must be a list")
    if (length(traps) != length(caphist))
        stop("number of sessions in traps does not match caphist")
    check.dim <- sapply(traps, nrow)
    if (!all(check.dim == caphist.dimensions[2, ]))
        stop("number of traps does not match caphist")
    if (!is.null(trapCovs)) {
        if (!is.list(trapCovs))
            stop("trapCovs must be a list")
        if (any(!sapply(trapCovs, is.list)))
            stop("trapCovs must be a list of lists")
        if (any(!unlist(sapply(trapCovs, function(x) sapply(x,
            is.data.frame)))))
            stop("trapCovs must be a list of dataframes")
        if (length(trapCovs) != length(caphist))
            stop("number of sessions in trapCovs does not match caphist")
        check.dim <- lapply(trapCovs, function(x) sapply(x, nrow))
        for (i in 1:length(check.dim)) {
            if (!all(check.dim[[i]] == caphist.dimensions[2,
                i]))
                stop("number of traps does not match caphist")
        }
    }
    if (!is.null(sigCovs)) {
        if (nrow(sigCovs) != length(caphist))
            stop("number of rows in sigCovs does not match number of sessions")
        if (!"session" %in% colnames(sigCovs)) {
            sigCovs$session <- factor(1:n.sessions)
        }
        if (!is.null(indCovs)) {
            if ("sex" %in% colnames(indCovs[[1]])) {
                sigCovs <- sigCovs[rep(1:n.sessions, 2), , drop = F]
                rownames(sigCovs) <- NULL
                sigCovs$sex <- factor(rep(c("female", "male"),
                  each = n.sessions))
            }
        }
    }
    else {
        sigCovs <- data.frame(session = factor(1:n.sessions))
        if (!is.null(indCovs)) {
            if ("sex" %in% colnames(indCovs[[1]])) {
                sigCovs <- sigCovs[rep(1:n.sessions, 2), , drop = F]
                rownames(sigCovs) <- NULL
                sigCovs <- sigCovs[rep(1:n.sessions, 2), , drop = F]
                rownames(sigCovs) <- NULL
                sigCovs$sex <- factor(rep(c("female", "male"),
                  each = n.sessions))
            }
        }
    }
    if (!is.null(trapOperation)) {
        if (!is.list(trapOperation))
            stop("trapOperation must be a list")
        if (length(trapOperation) != length(caphist))
            stop("number of sessions in trapOperation does not match caphist")
        check.dim <- sapply(trapOperation, nrow)
        if (!all(check.dim == caphist.dimensions[2, ]))
            stop("number of traps does not match caphist")
    }
    max.dist <- NULL
    for (i in 1:length(caphist)) {
        for (j in 1:nrow(caphist[[i]])) {
            if (dim(caphist[[i]])[3] > 1) {
                where <- apply(caphist[[i]][j, , ], 1, sum) >
                  0
            }
            else {
                where <- caphist[[i]][j, , ] > 0
            }
            if (sum(where) > 1)
                max.dist <- c(max.dist, max(0, dist(traps[[i]][where,
                  c("X", "Y")]), na.rm = T))
        }
    }
    mmdm <- mean(max.dist[max.dist > 0], na.rm = T)
    mdm <- max(max.dist, na.rm = T)
    if (!is.null(telemetry)) {
        if (!is.list(telemetry$fixfreq))
            stop("telemetry$fixfreq must be a list")
        fixfreq.dimensions <- sapply(telemetry$fixfreq, dim)
        if (nrow(fixfreq.dimensions) == 2)
            fixfreq.dimensions <- rbind(fixfreq.dimensions, 1)
        if (!is.null(telemetry$indCovs)) {
            if (!is.list(telemetry$indCovs))
                stop("telemetry$indCovs must be a list")
            if (any(!sapply(telemetry$indCovs, is.data.frame)))
                stop("telemetry$indCovs must be a list of dataframes")
            if (length(telemetry$indCovs) != length(telemetry$fixfreq))
                stop("number of sessions in telemetry$indCovs does not match telemetry$fixfreq")
            check.dim <- sapply(telemetry$indCovs, nrow)
            if (any(check.dim != fixfreq.dimensions[1, ]))
                stop("number of individuals in telemetry$indCovs does not match telemetry$fixfreq")
            if (any(!names(indCovs[[1]]) %in% c(names(telemetry$indCovs[[1]]),
                "removed")))
                stop("indCovs do not match between capture and telemetry data")
        }
        if (!is.null(telemetry$cap.tel)) {
            if (!is.list(telemetry$cap.tel))
                stop("telemetry$indCovs must be a list")
            warning("make sure captured individuals w/ collars sorted first!")
        }
    }
    if (!is.null(rsfDF)) {
        library(FNN)
        rsfCovs <- names(rsfDF[[1]][, -c(1:2), drop = F])
        if (is.null(trapCovs)) {
            trapCovs <- list()
            length(trapCovs) <- n.sessions
            for (s in 1:n.sessions) {
                trap.grid <- as.vector(get.knnx(rsfDF[[s]][,
                  c("X", "Y")], traps[[s]][, c("X",
                  "Y")], 1)$nn.index)
                trapCovs[[s]] <- list()
                length(trapCovs[[s]]) <- caphist.dimensions[3,
                  s]
                for (k in 1:caphist.dimensions[3, s]) {
                  trapCovs[[s]][[k]] <- data.frame(rsfDF[[s]][trap.grid,
                    rsfCovs])
                  names(trapCovs[[s]][[k]]) <- rsfCovs
                }
            }
        }
        else {
            for (s in 1:n.sessions) {
                if (any(!rsfCovs %in% trapCovs[[s]][[1]])) {
                  miss.rsfCovs <- rsfCovs[which(!rsfCovs %in%
                    trapCovs[[s]][[1]])]
                  trap.grid <- as.vector(get.knnx(rsfDF[[s]][,
                    c("X", "Y")], traps[[s]][, c("X",
                    "Y")], 1)$nn.index)
                  for (k in 1:caphist.dimensions[3, s]) {
                    newtrapCovs <- data.frame(rsfDF[[s]][trap.grid,
                      miss.rsfCovs])
                    names(newtrapCovs) <- miss.rsfCovs
                    trapCovs[[s]][[k]] <- data.frame(trapCovs[[s]][[k]],
                      newtrapCovs)
                  }
                }
            }
        }
    }
    scrFrame <- list(caphist = caphist, traps = traps, indCovs = indCovs,
        trapCovs = trapCovs, sigCovs = sigCovs, trapOperation = trapOperation,
        occasions = caphist.dimensions[3, ], type = type, mmdm = mmdm,
        mdm = mdm, telemetry = telemetry)
    class(scrFrame) <- "scrFrame"
    return(scrFrame)
}



##############################  SITE 4 #########################################

## MUESTREADOR PERSONALIZADO METROPOLIS-HASTINGS
## muestreador para actualizar conjuntamente y.un[1:M,j] de manera que ponemos
## en cada paso del muestreo la condicion de que sumen n[j]
IDSampler1 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # Definimos los componentes a usar
    nnidd1<-control$nnidd1
    j<-control$j
    M1<-control$M1
    calcNodes <- model$getDependencies(target)
  },

run = function() {
  lam.curr <- model$lamS1[1:M1,j] # Conteos esperados de individuos por trampa

  #Muestrea y[1:M,j] reasignando n[j] usando el condicional completo
  switch.probs <- lam.curr[1:M1]/sum(lam.curr[1:M1])

  #propone nuevas identificaciones para nnid[j,k]
  y.latent.curr <- model$y.fullS1[1:M1,j]- model$y.obs1[1:M1,j]
  y.latent.prop <- rmulti(1, nnidd1, switch.probs[1:M1])
  model$y.fullS1[1:M1,j] <<-  model$y.obs1[1:M1,j] + y.latent.prop

  # modelo inicial logProb
  model_lp_initial <- model$getLogProb(calcNodes)

  # modelo propuesto logProb
  model_lp_proposed <- model$calculate(calcNodes)

  # Relación log-Metropolis-Hastings
  log_MH_ratio <- (model_lp_proposed + dmulti(y.latent.curr, nnidd1, switch.probs, log=TRUE)) -
                  (model_lp_initial + dmulti(y.latent.prop,  nnidd1, switch.probs, log=TRUE))

  # Paso Metrópolis-Hastings
  accept <- decide(log_MH_ratio)
    if(accept) {
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    } else {
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    }
  },
  methods = list( reset = function () {} )
)

##############################  SITE 4 #########################################

## MUESTREADOR PERSONALIZADO METROPOLIS-HASTINGS
## muestreador para actualizar conjuntamente y.un[1:M,j] de manera que ponemos
## en cada paso del muestreo la condicion de que sumen n[j]
IDSampler2 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # Definimos los componentes a usar
    nnidd2<-control$nnidd2
    j<-control$j
    M2<-control$M2
    calcNodes <- model$getDependencies(target)
  },

run = function() {
  lam.curr <- model$lamS2[1:M2,j] # Conteos esperados de individuos por trampa

  #Muestrea y[1:M,j] reasignando n[j] usando el condicional completo
  switch.probs <- lam.curr[1:M2]/sum(lam.curr[1:M2])

  #propone nuevas identificaciones para nnid[j,k]
  y.latent.curr <- model$y.fullS2[1:M2,j]- model$y.obs2[1:M2,j]
  y.latent.prop <- rmulti(1, nnidd2, switch.probs[1:M2])
  model$y.fullS2[1:M2,j] <<-  model$y.obs2[1:M2,j] + y.latent.prop

  # modelo inicial logProb
  model_lp_initial <- model$getLogProb(calcNodes)

  # modelo propuesto logProb
  model_lp_proposed <- model$calculate(calcNodes)

  # Relación log-Metropolis-Hastings
  log_MH_ratio <- (model_lp_proposed + dmulti(y.latent.curr, nnidd2, switch.probs, log=TRUE)) -
                  (model_lp_initial + dmulti(y.latent.prop,  nnidd2, switch.probs, log=TRUE))

  # Paso Metrópolis-Hastings
  accept <- decide(log_MH_ratio)
    if(accept) {
      copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    } else {
      copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
    }
  },
  methods = list( reset = function () {} )
)

