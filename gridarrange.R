FIG3ABCD = grid.arrange(P2A, P2B, P2C, P2D, ncol = 2, nrow=2)

ggsave("Fig3ABCD.pdf",
       plot = FIG3ABCD,
       units="cm",
       width=15,
       height=15,
       dpi = 300)
