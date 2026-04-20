contibutionPlot.data <-
  function(hdp, clone.key) {
    contribution <- t(hdp$contr)
    absolute <- t(hdp$activities)

    tb <-
      Reduce(
        left_join,
        list(
          contribution %>%
            as.data.frame() %>%
            tibble::rownames_to_column("Signature") %>%
            tidyr::pivot_longer(-Signature,
              names_to = "Sample",
              values_to = "Contribution"
            ),
          t(hdp$contr_lb95) %>%
            as.data.frame() %>%
            tibble::rownames_to_column("Signature") %>%
            tidyr::pivot_longer(-Signature,
              names_to = "Sample",
              values_to = "Contribution_lb95"
            ),
          absolute %>%
            as.data.frame() %>%
            tibble::rownames_to_column("Signature") %>%
            tidyr::pivot_longer(-Signature,
              names_to = "Sample",
              values_to = "Absolute"
            )
        )
      )






    tb %<>% separate(
      col = Sample,
      into = c("Patient", "cluster_no"),
      convert = T,
      sep = "_",
      remove = F
    ) %>%
      relocate(Patient, .before = Sample) %>%
      relocate(Signature, .before = Contribution)


    # Annotate clone names & category
    # remove NA Clones
    tb <-
      left_join(
        tb,
        clone.key,
      ) %>%
      relocate(Clone_abv, .after = cluster_no) %>%
      relocate(category_label, .after = Patient)

    # add clone total
    tb <-
      left_join(
        tb,
        colSums(hdp$mutcount) %>%
          enframe(name = "Sample", value = "Sample_total")
      ) %>% relocate(Sample_total, .after = Sample)

    # Add csm

    csm <- hdp$sample_csm
    names(csm) <- colnames(hdp$mutcount)

    tb <-
      left_join(tb, csm %>% enframe(name = "Sample", value = "Sample_csm"))


    # Add Patient to Clone name and factor clones
    tb %<>% mutate(Clone_abv = paste0(Patient, "_", Clone_abv)) %>%
      dplyr::mutate(
        Sample = factor(Sample,
          levels = unique(Sample)
        ),
        Signature = factor(Signature,
          levels = unique(Signature)
        ),
        Clone_abv = factor(Clone_abv,
          levels = unique(Clone_abv)
        ),
        category_label = factor(
          category_label,
          levels = c(
            "Stable progressed",
            "Stable",
            "Transformed AML",
            "Transformed MF"
          )
        )
      )

    return(tb)
  }


contibutionPlot.con_plot <- function(tb, facet) {
  present_sigs <- tb %>%
    dplyr::filter(Contribution != 0) %>%
    dplyr::pull(Signature) %>%
    unique()

  tb %>% ggplot(aes(x = Clone_abv, y = Contribution, fill = Signature)) +
    geom_bar(
      aes(colour = Signature),
      position = "fill",
      stat = "identity",
      # colour = "black",
      linewidth = 0
    ) +
    labs(x = "", y = "Relative contribution") +
    scale_fill_discrete(breaks = present_sigs) +
    theme_bw() +
    theme(
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    coord_flip() +
    facet_grid(
      rows = vars(.data[[facet]]),
      scales = "free_y",
      space = "free"
    ) +
    theme(
      strip.background = element_blank(),
      strip.text.y = element_blank()
    ) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
}



contibutionPlot.abs_plot <- function(tb, facet) {
  present_sigs <- tb %>%
    dplyr::filter(Contribution != 0) %>%
    dplyr::pull(Signature) %>%
    unique()

  tb %>% ggplot(aes(x = Clone_abv, y = Absolute, fill = Signature)) +
    geom_bar(
      aes(colour = Signature),
      stat = "identity",
      # colour = "black",
      linewidth = 0
    ) +
    labs(x = "", y = "Absolute contribution \n (no. mutations)") +
    scale_fill_discrete(breaks = present_sigs) +
    theme_bw() +
    theme(
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    coord_flip() +
    facet_grid(
      rows = vars(.data[[facet]]),
      scales = "free_y",
      space = "free"
    ) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    scale_y_continuous(limits = c(0, max(pretty.default(
      range(tb %>% group_by(Clone_abv) %>% summarise(x = sum(Absolute)) %>% .$x)
    ))), expand = c(0, 0))
}
