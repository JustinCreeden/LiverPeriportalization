run_progeny = function (E, M, gene_name = "gene", value_name = "expression",
                        id_name = "sample", permutation = 10000, ...) {
  plan(multiprocess)
  E = E %>% mutate_if(is.factor, as.character)
  
  if (permutation > 0) {
    null_model = future_map_dfr(1:permutation, .progress = T, function(p) {
      E %>%
        group_by(!!!syms(id_name)) %>%
        sample_frac() %>%
        ungroup() %>%
        mutate(!!gene_name := E[[gene_name]]) %>%
        run_progeny(M, gene_name = gene_name, value_name = value_name,
                    id_name = id_name, permutation = 0)
    }) %>%
      group_by(!!!syms(id_name), pathway) %>%
      summarise(m = mean(activity),
                s = sd(activity)) %>%
      ungroup()
  }
  
  meta_data = E %>%
    select(-c(!!gene_name, !!value_name)) %>%
    distinct()
  
  emat = E %>%
    select(!!gene_name, !!id_name, !!value_name) %>%
    spread(!!id_name, !!value_name, fill = 0) %>%
    drop_na() %>%
    data.frame(row.names = 1, stringsAsFactors = F, check.names = F)
  
  model = M %>%
    spread(pathway, weight, fill = 0) %>%
    data.frame(row.names = 1, check.names = F, stringsAsFactors = F)
  
  common_genes = intersect(rownames(emat), rownames(model))
  emat_matched = emat[common_genes, , drop = FALSE] %>%
    t()
  model_matched = model[common_genes, , drop = FALSE] %>%
    data.matrix()
  
  stopifnot(names(emat_matched) == rownames(model_matched))
  
  progeny_scores = emat_matched %*% model_matched %>%
    data.frame(stringsAsFactors = F, check.names = F) %>%
    rownames_to_column(id_name) %>%
    gather(key = pathway, value = activity, -!!id_name) %>%
    as_tibble() %>%
    inner_join(meta_data, by = id_name)
  
  if (permutation > 0) {
    progeny_z_scores = progeny_scores %>%
      inner_join(null_model, by = c(id_name, "pathway")) %>%
      mutate(activity = (activity - m)/s) %>%
      select(!!id_name, pathway, activity)
    return(progeny_z_scores)
  }
  else {
    return(progeny_scores)
  }
}

run_viper = function(E, regulon, gene_name = "gene", value_name = "expression",
                     id_name = "sample", regulator_name = "tf",  ...) {
  meta_data = E %>%
    select(-c(!!gene_name, !!value_name)) %>%
    distinct()
  
  meta_regulon_data = regulon %>%
    select(-c(target, mor, likelihood)) %>%
    distinct()
  
  emat = E %>%
    select(!!gene_name, !!id_name, !!value_name) %>%
    spread(!!id_name, !!value_name, fill=0) %>%
    drop_na() %>%
    data.frame(row.names = 1, stringsAsFactors = F, check.names = F)
  
  viper_regulon = regulon %>%
    df2regulon(regulator_name = regulator_name)
  
  activity_scores = viper(eset = emat, regulon = viper_regulon, nes = T,
                          method = 'none', minsize = 4, eset.filter = F,
                          adaptive.size = F) %>%
    data.frame(stringsAsFactors = F, check.names = F) %>%
    rownames_to_column(var = regulator_name) %>%
    gather(key=!!id_name, value="activity", -!!regulator_name) %>%
    as_tibble() %>%
    inner_join(., meta_data, by=id_name) %>%
    inner_join(., meta_regulon_data, by = regulator_name)
  
  return(activity_scores)
}

df2regulon = function(df, regulator_name = "tf") {
  regulon = df %>%
    split(.[regulator_name]) %>%
    map(function(dat) {
      targets = setNames(dat$mor, dat$target)
      likelihood = dat$likelihood
      list(tfmode = targets, likelihood = likelihood)
    })
  return(regulon)
}