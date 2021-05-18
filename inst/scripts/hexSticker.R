devtools::load_all(".")
library("hexSticer")

imgurl <- system.file("images", "jacksparrow-01.png", package = "sparrow")
sticker(imgurl, package = "sparrow",
        p_size = 10,
        s_x = 1.05,
        s_y = 0.7,
        # s_height = 2,
        s_width = 0.55,
        p_color = "black",
        h_fill = "white",
        h_color = "black",
        filename = "man/figures/sparrow.png")
