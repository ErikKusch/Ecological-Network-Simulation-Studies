#' ####################################################################### #
#' PROJECT: [InfVal; Master Script for Parametrising Entire Analyses]
#' CONTENTS:
#'  - Master File for Framework Execution
#'  DEPENDENCIES:
#'  -
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
rm(list = ls())
set.seed(42)

# PREAMBLE ================================================================
## Directories ------------------------------------------------------------
message("Registering Directories")
Dir.Base <- file.path(getwd(), "Demo of NetDemoSims Capabilities")
Dir.Data <- file.path(Dir.Base, "Data")
Dir.Exports <- file.path(Dir.Base, "Exports")
Dirs <- c(Dir.Data, Dir.Exports)
CreateDir <- sapply(Dirs, function(x) if (!dir.exists(x)) dir.create(x))

## Packages ---------------------------------------------------------------
message("Loading Libraries")
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}
package_vec <- c(
  ## Result Packages
  "ggplot2", # for plotting
  "reshape2",
  "cowplot",
  "scales",
  "viridis",
  "ggnewscale",
  "ggpubr"
)
sapply(package_vec, install.load.package)

## NetDemoSims Setup ------------------------------------------------------
if ("NetDemoSims" %in% rownames(installed.packages()) == FALSE) {
  devtools::install_github("https://github.com/ErikKusch/NetDemoSims")
}
library(NetDemoSims)
package_vec <- c(package_vec, "NetDemoSims")

## Functionality -------------------------------------------------------
`%nin%` <- Negate(`%in%`)

label_row <- function(text) {
    ggplot() +
        geom_rect(
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "white",
            color = NA
        ) +
        annotate(
            "text",
            x = 0, y = 0,
            label = text,
            hjust = 0,
            fontface = "bold",
            size = 5
        ) +
        coord_cartesian(xlim = c(0, 1), ylim = c(-1, 1), clip = "off") +
        theme_void() +
        theme(
            plot.margin = margin(0, 0, 0, 0)
        )
}

# RUN SET-UP ===========================================================
## Simulation Parameters -----------------------------------------------
## global simulation parameters -----
n_runs <- 1e3 # number of networks to simulate and infer for

## Network Creation
n_spec <- 15 # number of species per network
NetworkType <- "Association" # type of links
Sparcity <- 0.1 # how many % of associations should be exactly 0
MaxStrength <- 1 # absolute maximum of interspecific links

## Initial Individual Creation
n_individuals <- 1e3 # number of individuals for initialisation
n_mode <- "each" # how to interpret the above number
Env_range <- c(0, 20) # environmental landscape range
Trait_sd <- 0.25 # standard deviation of traits per species

## Space Creation
Env_col <- Env_row <- 2e2
x_gradient <- function(x) x # function relating x-location to environmental value
y_gradient <- function(y) 0 # function relating y-location to environmental value

## Carrying Capacity Creation
k_range <- c(300, 300) # lower and upper limit of carrying capacities

## Simulation Parameters
d0 <- 0.35 # base death rate
b0 <- 0.65 # base birth rate
t_max <- 1e2 # simulation time
t_inter <- 1 # when to record data
Env_sd <- 1 # environmental maladaption SD, higher = more permissive environment
migration <- 0.5 # sd of 0 centred normal for relocation of offspring
migration_top <- 0.25 # point at which migration peaks
migration_trunc <- 0.75 # maximum distance of migration
Effect_Dis <- 0.5 # maximum distance at which link effect manifests
verbose <- TRUE # whether to produce simulation progress tracker in console

# NO NETWORK SPECIES IN ENVIRONMENT =========================================
Network_mat <- Sim.Network( # these are turned off in the actual simulation below
    n_spec = 4,
    NetworkType = "Association",
    Sparcity = 0,
    MaxStrength = 1, # what absolute value is the maximum link strength
    seed = 42 # seed for randomness
)
Plot.NetMat(Network_mat)

CarryingK_vec <- c(Sp_01 = 400, Sp_02 = 400, Sp_03 = 400, Sp_04 = 400)
Niches_vec <- c(Sp_01 = 2, Sp_02 = 4, Sp_03 = 6, Sp_04 = 8)

Initialise_df <- Sim.Initialise(
    n_spec = 4,
    n_individuals = 3e3,
    n_mode = "each",
    Env_range = c(0, 10),
    Trait_means = Niches_vec,
    Trait_sd = 0,
    seed = 42
)
Plot.IndivsInSpace(Initialise_df)

Env_mat <- Sim.Space(
    x_range = c(0, 10),
    y_range = c(0, 10),
    ncol = 1e3, nrow = 1e3,
    x_gradient = function(x) x,
    y_gradient = function(y) y
) / 2
Plot.Environment(Env_mat)

FNAME <- file.path(Dir.Data, "BoxDEMO_SpaceOnly_SimResult.RData")
if (file.exists(FNAME)) {
    load(FNAME)
} else {
    SimResult <- Sim.Compute(
        # Demographic parameters
        d0 = d0,
        b0 = b0,
        k_vec = CarryingK_vec,
        ID_df = Initialise_df,

        # Spatial parameters
        env.xy = Env_mat,
        env.sd = 0.5,
        mig.sd = 0.5,
        mig.top = 0.35,
        mig.trunc = 1,

        # Interaction parameters
        interac.maxdis = 0,
        interac.mat = Network_mat,
        interac.scale = 0,

        # Simulation parameters
        Sim.t.max = t_max * 2,
        Sim.t.inter = t_inter,
        seed = 0,
        verbose = verbose, # whether to print progress in time as current time
        RunName = "BoxDEMO_SpaceOnly",
        writeFile = FALSE
    )
    save(SimResult, file = FNAME)
}
last_df <- SimResult[[length(SimResult)]]

p_main <- Plot.Environment(Env_mat) +
    geom_point(
        data = last_df,
        aes(x = X, y = Y, shape = Species),
        fill = "white"
    ) +
    scale_shape_manual(values = 1:length(unique(last_df$Species)))

p_shape <- ggplot(
    last_df,
    aes(x = X, y = Y, shape = Species)
) +
    geom_point() +
    scale_shape_manual(values = 1:length(unique(last_df$Species))) +
    theme_void() +
    theme(legend.position = "right")

shape_legend <- get_legend(p_shape)

p_main_clean <- p_main +
    guides(shape = "none")

Enviro_SpecinSpace_gg <- plot_grid(
    p_main_clean,
    shape_legend,
    ncol = 2,
    rel_widths = c(1, 0.1)
)
ggsave(
    filename = file.path(Dir.Exports, "Figure_BoxDEMO_SpaceOnly_SpeciesInSpace.jpg"),
    plot = Enviro_SpecinSpace_gg,
    width = 18, height = 18, units = "cm"
)
print(Enviro_SpecinSpace_gg)

# 1-1 ASSOCIATIONS WITH ENVIRONMENT ========================================
AssocPanels <- list("Competition" = -1, "Mutualism" = 1, "No Association" = 0)

AssocPanels_gg <- lapply(1:length(AssocPanels), function(x) {
    x <- AssocPanels[[x]]
    if (x == 1) {
        y <- c("Normal", "OvercomingEnvironment")
    } else {
        y <- "Normal"
    }
    # print(x)

    plots <- lapply(y, function(i) {
        FNAME <- file.path(Dir.Data, paste0("BoxDEMO_AssocsPanels_SimResult_Assoc_", x, gsub(i, pattern = "Normal", replacement = ""), ".RData"))
        # print(FNAME)

        Network_mat <- Sim.Network(
            n_spec = 2,
            NetworkType = "Association",
            Sparcity = 0,
            MaxStrength = 1, # what absolute value is the maximum link strength
            seed = 42 # seed for randomness
        )
        Network_mat[2, 1] <- Network_mat[1, 2] <- x

        CarryingK_vec <- rep(200, ncol(Network_mat))
        names(CarryingK_vec) <- colnames(Network_mat)

        if (i == "OvercomingEnvironment") {
            Niches_vec <- c(Sp_01 = 0.5, Sp_02 = 1.5)
        } else {
            Niches_vec <- c(Sp_01 = 0, Sp_02 = 0)
        }

        Initialise_df <- Sim.Initialise(
            n_spec = 2,
            n_individuals = 2e3,
            n_mode = "each",
            Env_range = c(0, 10),
            Trait_means = Niches_vec,
            Trait_sd = 0,
            seed = 42
        )
        # Plot_IndivsInSpace(Initialise_df)

        if (i == "OvercomingEnvironment") {
            Env_mat <- Sim.Space(
                x_range = c(0, 10),
                y_range = c(0, 10),
                ncol = 1e3, nrow = 1e3,
                x_gradient = function(x) x,
                y_gradient = function(y) 0
            ) / 5
        } else {
            Env_mat <- Sim.Space(
                x_range = c(0, 10),
                y_range = c(0, 10),
                ncol = 1e3, nrow = 1e3,
                x_gradient = function(x, alpha = 1, L = 10) {
                    exp(-alpha * x) + exp(-alpha * (L - x))
                },
                y_gradient = function(y, alpha = 1, L = 10) {
                    exp(-alpha * y) + exp(-alpha * (L - y))
                }
            )
        }
        # Plot_Environment(Env_mat)

        if (file.exists(FNAME)) {
            load(FNAME)
        } else {
            SimResult <- Sim.Compute(
                # Demographic parameters
                d0 = d0,
                b0 = b0,
                k_vec = CarryingK_vec,
                ID_df = Initialise_df,

                # Spatial parameters
                env.xy = Env_mat,
                env.sd = ifelse(i == "OvercomingEnvironment", 0.5, 0.5),
                mig.sd = 0.5,
                mig.top = 0.25,
                mig.trunc = 1,

                # Interaction parameters
                interac.maxdis = 0.5,
                interac.mat = Network_mat,
                interac.scale = 1,

                # Simulation parameters
                Sim.t.max = t_max,
                Sim.t.inter = t_inter,
                seed = 42,
                verbose = verbose, # whether to print progress in time as current time
                RunName = "BoxDEMO_AssocsPanels",
                writeFile = FALSE
            )
            save(SimResult, file = FNAME)
        }

        last_df <- SimResult[[length(SimResult)]]

        p_main <- Plot.Environment(Env_mat) +
            geom_point(
                data = last_df,
                aes(x = X, y = Y, shape = Species),
                fill = "white", col = "white"
            ) +
            scale_shape_manual(values = c(3, 5)) +
            theme(legend.position = "top")
        p_main
    })
    plots
})
names(AssocPanels_gg) <- names(AssocPanels)

legend <- get_legend(AssocPanels_gg[[1]][[1]])

AssocPanels_gg <- unlist(AssocPanels_gg[c(3, 1, 2)])

AssocPanels_gg2 <- plot_grid(
    legend,
    plot_grid(plotlist = lapply(AssocPanels_gg, FUN = function(x) {
        x + theme(legend.position = "none") + labs(x = "", y = "")
    }), ncol = 2),
    ncol = 1,
    rel_heights = c(0.1, 1)
)
ggsave(
    filename = file.path(Dir.Exports, "Figure_BoxDEMO_AssocPanels_gg.jpg"),
    plot = AssocPanels_gg2,
    width = 14 * 1.7, height = 16 * 1.5, units = "cm"
)

# ALL-IN-ONE ASSOCIATIONS WITH ENVIRONMENT =================================
FNAME <- file.path(Dir.Data, "BoxDEMO_AssocsAllinOne_SimResult.RData")

Network_mat <- Sim.Network(
    n_spec = 6,
    NetworkType = "Association",
    Sparcity = 0,
    MaxStrength = 0, # what absolute value is the maximum link strength
    seed = 42 # seed for randomness
)
Network_mat[2, 1] <- Network_mat[1, 2] <- 0.5
Network_mat[5, 6] <- Network_mat[6, 5] <- -0.5

CarryingK_vec <- rep(300, ncol(Network_mat))
names(CarryingK_vec) <- colnames(Network_mat)
Niches_vec <- c(Sp_01 = 5, Sp_02 = 5, Sp_03 = 10, Sp_04 = 10, Sp_05 = 15, Sp_06 = 15)

Initialise_df <- Sim.Initialise(
    n_spec = 6,
    n_individuals = 2e3,
    n_mode = "each",
    Env_range = c(0, 20),
    Trait_means = Niches_vec,
    Trait_sd = 0,
    seed = 42
)
Plot.IndivsInSpace(Initialise_df)

Env_mat <- Sim.Space(
    x_range = c(0, 20),
    y_range = c(0, 20),
    ncol = 1e3, nrow = 1e3,
    x_gradient = function(x) x,
    y_gradient = function(y) 0
)
Plot.Environment(Env_mat)

if (file.exists(FNAME)) {
    load(FNAME)
} else {
    SimResult <- Sim.Compute(
        # Demographic parameters
        d0 = d0,
        b0 = b0,
        k_vec = CarryingK_vec,
        ID_df = Initialise_df,

        # Spatial parameters
        env.xy = Env_mat,
        env.sd = 1,
        mig.sd = 0.75,
        mig.top = 0.5,
        mig.trunc = 1,

        # Interaction parameters
        interac.maxdis = 1,
        interac.mat = Network_mat,
        interac.scale = 1,

        # Simulation parameters
        Sim.t.max = t_max,
        Sim.t.inter = t_inter,
        seed = 42,
        verbose = verbose, # whether to print progress in time as current time
        RunName = "BoxDEMO_Assocs",
        writeFile = FALSE
    )
    save(SimResult, file = FNAME)
}
last_df <- SimResult[[length(SimResult)]]

p_main <- Plot.Environment(Env_mat) +
    geom_point(
        data = last_df,
        aes(x = X, y = Y, shape = Species),
        fill = "white", col = "white"
    ) +
    scale_shape_manual(values = 1:length(unique(last_df$Species)))

p_shape <- ggplot(
    last_df,
    aes(x = X, y = Y, shape = Species)
) +
    geom_point() +
    scale_shape_manual(values = 1:length(unique(last_df$Species))) +
    theme_void() +
    theme(legend.position = "right")

shape_legend <- get_legend(p_shape)

p_main_clean <- p_main +
    guides(shape = "none")

SpecinSpace2_gg <- plot_grid(
    p_main_clean,
    shape_legend,
    ncol = 2,
    rel_widths = c(1, 0.1)
)

ggsave(
    filename = file.path(Dir.Exports, "Figure_BoxDEMO_Assoc_AllinOne_gg.jpg"),
    plot = SpecinSpace2_gg,
    width = 14 * 1.7, height = 16 * 1.5, units = "cm"
)


# INTERACTIONS WITH ENVIRONMENT ============================================
Env_mat <- Sim.Space(
    x_range = c(0, 10),
    y_range = c(0, 10),
    ncol = 1e3, nrow = 1e3,
    x_gradient = function(x) {
        2 * ((x - 5)^2) / 2
    },
    y_gradient = function(y) {
        2 * ((y - 5)^2) / 2
    }
) / 5
Env_mat <- abs(Env_mat - 3)
Plot.Environment(Env_mat)

stop("make trophic network")
# maybe make one unaffected primary, one secondary consumer, and one tertiary predator? maybe with weaker links, specify ncihe preference of primary and secondary consumer as increasingly unaligned with environment (do this via negative preferences)
