process_wgcna_module<- function(dir){
  path <- list.files(path = dir, pattern = "wgcna_module.rds", full.names = T)
  print(path)
  wgcna_module <- readRDS(path)
  
  module_list <- list("original_wgcna_module" = wgcna_module)
  
  unadjusted_p_val_module <- MODifieRDev::wgcna_adjust_significance(p_value = 0.05, wgcna_module = wgcna_module, use_unadjusted = T)
  
  module_list <- c(module_list, "unadjusted_p_value" = list(unadjusted_p_val_module))
  if (length(wgcna_module$module_colors != 0)){
    tryCatch(modules_by_color <- MODifieRDev::wgcna_split_module_by_color(wgcna_module = wgcna_module))

      names(modules_by_color) <- wgcna_module$module_colors
    
    module_list <- c(module_list, modules_by_color)
    
  }
  steps <- seq(from = 200, to = 1200, by = 200)
  
  modules_by_size <- lapply(X = steps, FUN = module_by_size, wgcna_module = wgcna_module, output_file = paste0(dir, "/wgcna_size_module_composition.csv"))
  names(modules_by_size) <- paste0("size_", as.character(steps))
  
  module_list <- c(module_list, modules_by_size)
  
  write_info_table(wgcna_module = wgcna_module, output_file = paste0(dir, "/wgcna_information_table.csv"))
  write_module_list(module_list = module_list, output_file = paste0(dir, "/wgcna_module_list.csv"))
  
}

write_module_list <- function(module_list, output_file){
  for (i in 1:length(module_list)){
    invisible(write.table(x = t(c(names(module_list)[[i]], "", module_list[[i]]$module_genes)),
                          file = output_file, append = T, row.names = F,
                          col.names = F, sep = "\t", quote = F))
  }
}

write_info_table <- function(wgcna_module, output_file){
  info_table <- cbind(wgcna_module$correlation_to_trait_table, lengths(wgcna_get_all_module_genes(wgcna_module)))
  
  colnames(info_table)[4] <- "module_size"
  
  write.table(x = info_table, file = output_file, row.names = T, col.names = T)
}

module_by_size <- function(size, wgcna_module, output_file){
  size_module <- MODifieRDev::wgcna_set_module_size(size = size, wgcna_module = wgcna_module)
  write.table(x = t(c(as.character(size), size_module$module_colors)),row.names = F, col.names = F, file = output_file, append = T)
  return(size_module)
}


dirs <- list.dirs(path = "/home/dirk/wgcna_modules/v_10_5_update/")

dirs <- dirs[-1]

lapply(X = dirs, FUN = process_wgcna_module)

