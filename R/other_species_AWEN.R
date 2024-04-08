


#' Calculates the AWEN proportions of different component based on the data in Repo et al., 2014. The function is a collection of factors and can be used for both coarse and fine components by
#' selecting the specific option
#'
#' @param biom biomass to be converted into AWEN pool, unit does not matter (output will be in the same unit)
#' @param comp components. "coarse" = stem, branch, roots, stump; "fine" = foliage, fine roots
#' @param type tree type, "temp.broad.ev" = temperate broadleaved evergreen, "temp.broad.sum" =  = temperate broadleaved evergreen, bor.broad.sum = boreal broadleaved summergreen, "temp" = temperate, "bor.con" = boreal coniferous
#' @returns the mass of the four AWEN components
#' @examples
#' @references Repo, A., Böttcher, H., Kindermann, G., & Liski, J. (2014). sustainability of forest bioenergy in europe: land‐use‐related carbon dioxide emissions of forest harvest residues. GCB bioenergy, 7(4), 877-887. https://doi.org/10.1111/gcbb.12179
#' @seealso \link{foliage.AWEN, fineroot.AWEN, branches.AWEN}
#' @author Lorenzo Menichetti
#'
multispecies.repo.AWEN <- function(biom, comp, type) {

  biom.AWEN <- matrix(0,nrow=length(biom), ncol=4)

  # data frame with the awen factors
  repo_data <- data.frame(
    type = c("temp.broad.ev", "temp.broad.ev", "temp.broad.sum", "temp.broad.sum", "bor.broad.sum", "bor.broad.sum", "temp.con", "temp.con", "bor.con", "bor.con"),
    comp = c("coarse", "fine", "coarse", "fine", "coarse", "fine", "coarse", "fine", "coarse", "fine"),
    A = c(0.76, 0.49, 0.76, 0.39, 0.76, 0.39, 0.68, 0.51, 0.68, 0.5),
    W = c(0.01, 0.15, 0.01, 0.09, 0.01, 0.09, 0.02, 0.13, 0.01, 0.09),
    E = c(0.01, 0.08, 0.01, 0.05, 0.01, 0.05, 0.01, 0.10, 0.05, 0.05),
    N = c(0.22, 0.29, 0.22, 0.47, 0.22, 0.47, 0.29, 0.26, 0.26, 0.36),
    stringsAsFactors = FALSE
  )



  biom.AWEN[,1] <- repo_data[repo_data$type == type & repo_data$comp == comp,1]*biom

  biom.AWEN[,2] <- repo_data[repo_data$type == type & repo_data$comp == comp,1]*biom

  biom.AWEN[,3] <-repo_data[repo_data$type == type & repo_data$comp == comp,1]*biom

  biom.AWEN[,4] <- repo_data[repo_data$type == type & repo_data$comp == comp,1]*biom

  return(biom.AWEN)
}

