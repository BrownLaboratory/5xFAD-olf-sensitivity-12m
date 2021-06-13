library(cowplot)
library(here)


plot_grid(plotAcc12m 
          + geom_hline(yintercept = .85,
                       linetype = 'dashed')
          + coord_cartesian(ylim = c(.35, 1))
          + theme(legend.position = c(.8, .2)),
          plotAcc6m 
          + geom_hline(yintercept = .85,
                       linetype = 'dashed')
          + coord_cartesian(ylim = c(.35, 1))
          + theme(legend.position = 'none'),
          plotBestAcc12m, 
          plotBestAcc6m
          + theme(legend.position = 'none'),
          plot85Pr12m
          + theme(legend.position = 'none'), 
          plot85Pr6m
          + theme(legend.position = 'none'),
          nrow = 3,
          align = 'h',
          labels = 'AUTO'
)
ggsave(here('plots', 'newNewPlots', 'acc.png'),
       width = 8.5,
       height = 11)


plot_grid(plotDeePrimBox12m + coord_cartesian(ylim = c(0, 5)), 
          plotDeePrimBox6m + coord_cartesian(ylim = c(0, 5)),
          plotDeePrime12m + theme(legend.position = c(.8, .2)) + coord_cartesian(ylim = c(0, 5)), 
          plotDeePrime6m + theme(legend.position = 'none') + coord_cartesian(ylim = c(0, 5)),
          align = 'h',
          labels = 'AUTO')
ggsave(here('plots', 'newNewPlots', 'deePrime.png'),
       width = 8.5,
       height = 8.5)


plot_grid(plotResponseBiasBox12m + coord_cartesian(ylim = c(-1.5, .65)), 
          plotResponseBiasBox6m + coord_cartesian(ylim = c(-1.5, .65)),
          plotResponseBiasConcBox12m 
          + coord_cartesian(ylim = c(-1.8, 1))
          + theme(legend.position = c(.8, .2),
                  legend.background = element_blank()), 
          plotResponseBiasConcBox6m 
          + coord_cartesian(ylim = c(-1.8, 1))
          + theme(legend.position = 'none'), # c(.8, .2)),
          align = 'h',
          labels = 'AUTO')
ggsave(here('plots', 'newNewPlots', 'responseBias.png'),
       width = 8.5,
       height = 8.5)


# plot_grid(plotBestAcc12m, 
#           plotBestAcc6m
#           + theme(legend.position = 'none'),
#           plot85Pr12m
#           + theme(legend.position = 'none'), 
#           plot85Pr6m
#           + theme(legend.position = 'none'),
#           align = 'h',
#           labels = 'AUTO')
# ggsave(here('plots', 'newNewPlots', 'threshold.png'),
#        width = 8.5,
#        height = 8.5)

plot_grid(plotFa12m
          + theme(legend.position = c(.8, .85)), 
          plotFa6m
          + theme(legend.position = 'none'),
          align = 'h',
          labels = 'AUTO')
ggsave(here('plots', 'newNewPlots', 'faRate.png'),
       width = 8.5,
       height = 5.5)
