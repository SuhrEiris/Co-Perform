FIG3ABCD = grid.arrange(P2A, P2B, P2C, P2D, ncol = 2, nrow=2)

ggsave("Fig3ABCD.pdf",
       plot = FIG3ABCD,
       units="cm",
       width=15,
       height=15,
       dpi = 300)

FIGS2ABCD = grid.arrange(PS4A, PS4B, PS4C, noave_2, ncol = 2, nrow=2)

ggsave("FigS2ABCD.pdf",
       plot = FIGS2ABCD,
       units="cm",
       width=18,
       height=15,
       dpi = 300)
