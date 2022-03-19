readBlosum62<- function(){
  
  blosum62_txt <- file("Blosum62.txt",open="r")
  blosum62_table <-read.table(blosum62_txt,header=TRUE,check.names=FALSE)
  
  close(blosum62_txt)
  
  return (blosum62_table)
}

getBlosumScore <- function(Protein1, Protein2) {
  
  if(!any(grepl(Protein1, colnames(blosum62)))){
    Protein1<-"*"
  }
  
  if(!any(grepl(Protein2, colnames(blosum62)))){
    Protein2<-"*"
  }
  
  return(blosum62[Protein1,Protein2])
}

gapPenalty<-function(){
  
  options(max.print=999999)
  
  gap_penalty = readline(prompt = "Please enter gap penalty for alignments:");
  gap_penalty = as.integer(gap_penalty)
  
  return (gap_penalty)
}

readInput <- function(){
  
  fileName <- "Input3.txt"
  
  con <- file(fileName,open="r")
  
  seq_names <- c()
  seq_seq <- c()
  
  lin <-readLines(con)
  
  for (i in 1:length(lin)){
    if (i %% 2 == 0 ){
      lin[i] <- trimws(lin[i],"b")
      seq_seq <- append(seq_seq , lin[i])
    }
    else {
      lin[i] <- trimws(lin[i],"b")
      tmp <- lin[i]
      tmp <- gsub("[\r\n\t]", "", tmp)
      seq_names <- append(seq_names , tmp)
    }
  }
  
  close(con)
  
  inputs <- list()
  inputs <- append(inputs,list(seq_seq))
  inputs <- append(inputs,list(seq_names))
  
  return(inputs)
}

pairwiseAlignment <- function(seq1, seq2) {
  
  seq1 <- paste(" ", seq1 ,sep = "")
  seq2 <- paste(" ", seq2 , sep="" )
  
  
  
  arrow_table <- matrix(, nrow = nchar(seq1)  , ncol = nchar(seq2)  )
  score_table <- matrix(, nrow = nchar(seq1)  , ncol = nchar(seq2)  )
  
  
  seq1_chars <- strsplit(seq1, "")[[1]]
  seq2_chars <- strsplit(seq2, "")[[1]]
  
  
  colnames(arrow_table) <- seq2_chars
  colnames(score_table) <- seq2_chars
  
  rownames(arrow_table) <- seq1_chars
  rownames(score_table) <- seq1_chars
  
  for (i in 1:nrow(arrow_table)){
    for (j in 1:ncol(arrow_table)){
      if (i == 1 && j == 1 ){
        arrow_table[i,j] <- "."
        score_table[i,j] <- 0
      }
      else if (i == 1 && j != 1 ){
        arrow_table[i,j] = "<-"
        score_table[i,j] <- score_table[i,j-1] - gap_penalty 
      }
      else if (j == 1 && i != 1 ){
        arrow_table[i,j] = "^"
        score_table[i,j] <- score_table[i-1,j] - gap_penalty 
      }
      else {
        
        score <- (getBlosumScore(rownames(score_table)[i],colnames(score_table)[j]))
        
        if ( ( (score_table[i-1,j-1] + score  ) >= ( score_table[i,j-1] - gap_penalty ) ) && ((score_table[i-1,j-1] + score) >= (score_table[i-1,j] - gap_penalty ) ) ) {
          
          score_table[i,j] <- ( score_table[i-1,j-1] + score ) 
          arrow_table[i,j] <- "/"
          
        }
        
        else if ( ( (score_table[i-1,j] - gap_penalty ) >= ( score_table[i-1,j-1] + score ) ) && (score_table[i-1,j] -gap_penalty >= (score_table[i,j-1] - gap_penalty ) ) ) {
          
          score_table[i,j] <- ( score_table[i-1,j] - gap_penalty )
          arrow_table[i,j] <- "^"
          
        }
        else if (  ( (score_table[i,j-1] - gap_penalty ) >= ( score_table[i-1,j-1] + score) ) && (score_table[i,j-1] -gap_penalty >= (score_table[i-1,j] -gap_penalty ) ) ){
          
          score_table[i,j] <- ( score_table[i,j-1] - gap_penalty  )
          arrow_table[i,j] <- "<-"
          
        }
        
      }
      
    }
  }
  
  return(arrow_table)
  
}

backTracking <- function(arrow_table){
  
  col_seqs <- list()
  row_seqs <- list()
  result_seqs <- ""
  
  col_alignments <- strsplit(colnames(arrow_table), split = "")
  row_alignments <- strsplit(rownames(arrow_table), split = "")
  
  for ( i in 1:(lengths(col_alignments[2]))){
    col_seqs <- append(col_seqs,"")
    
  }
  
  for ( i in 1:(lengths(row_alignments[2]))){
    row_seqs <- append(row_seqs,"")
    
  }
  
  row = nrow(arrow_table)
  column= ncol(arrow_table)
  
  
  while(row!=1 || column !=1){
    
    if(any(grepl("/", arrow_table[row,column]))){
      
      for ( i in 1:(lengths(col_alignments[2]))){
        col_seqs[i] <- paste(col_alignments[[column]][i],col_seqs[i],sep="")
        
      }
      for ( i in 1:(lengths(row_alignments[2]))){
        row_seqs[i] <- paste(row_alignments[[row]][i],row_seqs[i],sep="")
        
      }
      
      row=row-1
      column=column-1
      
    }else if(any(grepl("<-", arrow_table[row,column]))){
      
      for ( i in 1:(lengths(col_alignments[2]))){
        col_seqs[i] <- paste(col_alignments[[column]][i],col_seqs[i],sep="")
        
      }
      for ( i in 1:(lengths(row_alignments[2]))){
        row_seqs[i] <- paste("_",row_seqs[i],sep="")
        
      }
      
      column=column-1
      
    }
    else if(any(grepl("^", arrow_table[row,column]))){
      
      for ( i in 1:(lengths(col_alignments[2]))){
        col_seqs[i] <- paste("_",col_seqs[i],sep="")
        
      }
      for ( i in 1:(lengths(row_alignments[2]))){
        row_seqs[i] <- paste(row_alignments[[row]][i],row_seqs[i],sep="")
        
      }
      
      row=row-1
    }
  }
  
  flag<-TRUE
  
  for(i in 1:(length(row_seqs))){
    
    if(flag){
      result_seqs <- paste(result_seqs,row_seqs[i],sep="")
      flag<-FALSE
    }
    else{
      result_seqs <- paste(result_seqs,row_seqs[i],sep="-")
    }
  }
  
  for(i in 1:(length(col_seqs))){
    result_seqs <- paste(result_seqs,col_seqs[i],sep="-")
  }
  return (result_seqs)
  
}

getSimilarityScore <- function(seq){
  
  sequences_tmp <- strsplit(seq, split = "-")
  sequences <- sequences_tmp [[1]] 
  
  s1 <- strsplit(sequences[1],split="")[[1]]
  
  s2 <- strsplit(sequences[2],split="")[[1]]
  
  matches <- 0 
  for ( i in 1: length(s1)){
    
    if (s1[i] == s2[i]){
      matches <- matches + 1
    }
    
  }
  return ( (matches / length(s1) * 100) ) 
}

allPairwiseAlignments <- function(seq_seq_list , seq_name_list){
  
  distance_matrix <- array(numeric(),c(length(seq_name_list),length(seq_name_list),1)) 
  
  colnames(distance_matrix) <- seq_name_list
  rownames(distance_matrix) <- seq_name_list
  
  for ( i in 1:length(colnames(distance_matrix))){
    for (j in 1:length(rownames(distance_matrix))){
      
      arw_tbl <- pairwiseAlignment(seq_seq_list[i] , seq_seq_list[j] )
      str <- backTracking(arw_tbl)
      val <- getSimilarityScore(str)
      distance_matrix[i,j,1] <- val
      
    }
  }
  return (distance_matrix)
}

buildGuideTree <- function(score_matrix){
  
  colnames <- colnames(score_matrix)
  rownames <- rownames(score_matrix)
  score_matrix <- matrix(score_matrix, nrow = nrow(score_matrix), ncol = ncol(score_matrix))
  
  for (i in 1:nrow(score_matrix)) {
    for (j in 1:ncol(score_matrix)){
      score_matrix[i, j] <- 100 - score_matrix[i, j] 
    }
  }
  
  colnames(score_matrix) <- colnames
  rownames(score_matrix) <- rownames
  
  
  fill_r_matrix <- function(){
    r_matrix <- matrix(
      nrow = 1,
      ncol = ncol(score_matrix)
    )
    
    sum <- 0
    
    for(row in 1:nrow(score_matrix)){
      for(col in 1:ncol(score_matrix)){
        if(row != col)
          sum = sum + score_matrix[row, col]
      }
      
      r_matrix[row] <- sum
      sum = 0
    }
    return(r_matrix)
  }
  
  
  fill_distance_matrix <- function(){
    
    distance_matrix <- matrix(
      nrow = nrow(score_matrix),
      ncol = ncol(score_matrix)
    )
    
    n <- 1
    
    for(row in 1:nrow(distance_matrix)){
      for(col in 1:n){
        if(row != col){
          distance_matrix[row, col] <- ((nrow(score_matrix) - 2) * score_matrix[row, col]) - (r_matrix[col] + r_matrix[row]) 
        }
      }
      n <- n + 1
    }
    return(distance_matrix)
  }
  
  
  choose_smallest <- function(){
    n <- 1
    min <- .Machine$integer.max
    rowMin <- 0
    colMin <- 0
    for(row in 1:nrow(distance_matrix)) {
      for(col in 1:n) {
        if(row != col && distance_matrix[row, col] < min){
          min <- distance_matrix[row, col]
          rowMin <- row
          colMin <- col
        }
      }
      n <- n + 1
    }
    return(c(colMin, rowMin))
  }
  
  
  edge_length <- function(minIndexes){
    colMin <- minIndexes[1]
    rowMin <- minIndexes[2]
    ColumnDis <- (score_matrix[rowMin, colMin] / 2) + ((r_matrix[colMin] - r_matrix[rowMin]) / (2 * (nrow(score_matrix) - 2)))
    RowDis <- score_matrix[rowMin, colMin] - ColumnDis
    return(paste(ColumnDis, RowDis, sep = "-") )
  }
  
  
  update_score_matrix <- function(){
    colMin <- minIndexes[1]
    rowMin <- minIndexes[2]
    new_scores <- function(value){
      return((score_matrix[colMin, value] + score_matrix[rowMin, value] - score_matrix[rowMin, colMin]) / 2)
    }
    
    new_score_matrix <- matrix(
      nrow = 1,
      ncol = ncol(score_matrix)
    )
    
    for(col in 1:ncol(new_score_matrix)){
      if(col != rowMin && col != colMin){
        new_score_matrix[col] <- new_scores(col)
      }
    }
    
    rowMinId <- rownames(score_matrix)[rowMin]
    colMinId <- colnames(score_matrix)[colMin]
    score_matrix[colMin, ] <- new_score_matrix
    score_matrix[, colMin] <- new_score_matrix
    score_matrix <- score_matrix[, -rowMin]
    score_matrix <- score_matrix[-rowMin, ]
    name_str <- paste(colMinId, rowMinId, sep = "-")
    rownames(score_matrix)[colMin] <- name_str
    colnames(score_matrix)[colMin] <- name_str
    score_matrix[colMin, colMin] <- 0
    
    return(score_matrix)
  }
  
  
  colnames_score_matrix <- colnames(score_matrix)
  tree <- list()
  
  while(dim(score_matrix) != c(2, 2)){
    
    r_matrix <- fill_r_matrix()
    distance_matrix <- fill_distance_matrix()
    minIndexes <- choose_smallest()
    index_str <- paste(colnames(score_matrix)[minIndexes[1]], colnames(score_matrix)[minIndexes[2]], sep = "-")
    length_str <- edge_length(minIndexes)
    tree <- append(tree, paste(index_str, length_str, sep = "/"))
    score_matrix <- update_score_matrix()
    colnames_score_matrix[minIndexes[1]] <- sprintf("%s%s%s%s%s", "(", colnames_score_matrix[minIndexes[1]], "-", colnames_score_matrix[minIndexes[2]], ")")
    colnames_score_matrix <- colnames_score_matrix[-minIndexes[2]]
    
  }
  
  cat("\n")
  print(score_matrix)
  cat("\n")
  index_str <- paste(colnames(score_matrix)[1], colnames(score_matrix)[2], sep = "-")
  length_str <- score_matrix[1, 2]
  tree <- append(tree, paste(index_str, length_str, sep = "/"))
  print(sprintf("%s%s%s%s", "Constructed Tree: ", colnames_score_matrix[1], "-", colnames_score_matrix[2]))
  print(tree)
  
  return (tree)
}

spScore <- function(vector_sp_1, vector_sp_2){
  
  vector_sp <- c(vector_sp_1, vector_sp_2) 
  
  sp_score<-0
  
  for ( i in 1:(length(vector_sp)-1)){
    for( j in (i+1):length(vector_sp)){
      if(vector_sp[i] == "_" && vector_sp[j]== "_"){
        sp_score<-sp_score+getBlosumScore(vector_sp[i],vector_sp[j])
      }
      else if(vector_sp[i] == "_" || vector_sp[j]== "_"){
        sp_score<-sp_score-gap_penalty
      }
      else{
        sp_score<-sp_score+getBlosumScore(vector_sp[i],vector_sp[j])
      }
    }
  }
  sp_score <- ((-1) * sp_score)
  return (sp_score)
}

progressiveAlignment <- function(gd_tree, seq_name_list , seq_fasta_list){
  
  alignments <- list()
  
  for (i in 1:length(gd_tree)){
    
    step <- strsplit(gd_tree[[i]], split = "/")[[1]]
    names <- strsplit(step[1],split = "-")[[1]]
    distances <- strsplit(step[2], split = "-")[[1]]
    temp <- ""
    
    if (length(names) == 2 ){
      seq1 <- ""
      seq2 <- ""
      flag <- 0 
      
      for (j in 1:length(names)){
        for (k in 1:length(seq_name_list)){
          if (flag == 0 && seq_name_list[k] == names[j] ) {
            flag <- 1 
            seq1 <- seq_fasta_list[k]
            break 
          }
          else if (flag == 1 && names[j] == seq_name_list[k]){
            seq2 <- seq_fasta_list[k]
            break 
          }
        }
      }
      
      temp <- step[1]
      alignment <- backTracking(pairwiseAlignment(seq1, seq2))
      temp <- paste(temp, alignment, sep="/")
      alignments <- append(alignments, temp)
      
    }
    else {
      seq1 <- ""
      seq2 <- ""
      flag_1 <- 0 
      flag_2 <- 0 
      index_flag <- 0
      index <- 1
      
      for (j in 1:length(names)){
        temp_index <- 0 
        
        for (k in 2: ( length(names) - 1 )) {
          temp_str <- ""
          
          for (t in 1:k ){
            if ( (j + (t-1) ) <= length(names) ) {
              if ( nchar (temp_str) == 0 ) {
                temp_str <- paste (temp_str , names[j+(t-1)],sep = "")
              }
              else {
                temp_str <- paste (temp_str , names[j+(t-1)],sep = "-")
              }
              
              temp_index <- (j + t)
              
            }
          }
          for (kk in  1:length(alignments)){
            cmp_tmp_step <- strsplit(alignments[[kk]], split = "/")[[1]]
            cmp_tmp <- cmp_tmp_step[1]
            
            if (cmp_tmp == temp_str && flag_1 == 0){
              alignments[[kk]] <- " "
              flag_1 <- 1
              seq1 <- cmp_tmp_step[2]
              index <- temp_index
              
              if (j != 1 ){
                index <- 1
              }
              break
            }
            else if (cmp_tmp == temp_str && flag_1 == 1){
              alignments[[kk]] <- " "
              flag_2 <- 1
              seq2 <- cmp_tmp_step[2]
              break
            }
            
          }
          
        }
        
      }
      if (flag_2 == 0 ){
        for (q in index:length(names)){
          for (p in 1:length(seq_name_list)) {
            if (names[q] == seq_name_list[p]){
              seq2 <- seq_fasta_list[p]
              if (index == 1 ){
                tmp <- seq2
                seq2 <- seq1
                seq1 <- tmp
              }
              break
            }
          }
          if (names[q] == seq_name_list[p]){
            break
          }
          
        }
      }
      
      #-----------------------------------------------
      seq1_sequences <- strsplit(seq1, split="-")[[1]]
      seq2_sequences <- strsplit(seq2, split="-")[[1]]

      seq1_to_align <- seq1_sequences[1] 
      seq2_to_align <- seq2_sequences[1] 
      
      arrow_table <- matrix(, nrow = nchar(seq1_to_align) + 1   , ncol = nchar(seq2_to_align)  + 1  )
      score_table <- matrix(, nrow = nchar(seq1_to_align) + 1   , ncol = nchar(seq2_to_align)  + 1   )
      
      seq1_chars <-  c("  ")
      seq2_chars <-  c("  ")
      
      for (a in 1:length(seq1_sequences)){
        tmp <- strsplit(seq1_sequences[a], split = "")[[1]]
        
        for (b in 1:length(tmp)){
          if (a != 1 ){
            seq1_chars[b+1] <- paste(seq1_chars[b+1], tmp[b], sep="")
          }
          else {
            seq1_chars <- append(seq1_chars, tmp[b])
          }
        }
      }
      
      for (a in 1:length(seq2_sequences)){
        tmp <- strsplit(seq2_sequences[a] , split = "")[[1]]
        
        for (b in 1:length(tmp)){
          if (a != 1 ){
            seq2_chars[b+1] <- paste(seq2_chars[b+1], tmp[b], sep="")
          }
          else {
            seq2_chars <- append(seq2_chars, tmp[b])
          }
        }
      }
      
      colnames(arrow_table) <- seq2_chars
      colnames(score_table) <- seq2_chars
      
      rownames(arrow_table) <- seq1_chars
      rownames(score_table) <- seq1_chars
      
      for (z in 1:nrow(arrow_table)){
        for (x in 1:ncol(arrow_table)){
          vector_sp_1_tmp <- strsplit(rownames(score_table)[z], split="")[[1]]
          vector_sp_2_tmp <- strsplit(colnames(score_table)[x], split ="") [[1]]
          
          if (z == 1 && x == 1 ){
            arrow_table[z,x] <- "."
            score_table[z,x] <- 0
          }
          else if (z == 1 && x != 1 ){
            arrow_table[z,x] = "<-"
            vector_sp_2 <- c()
            
            for (a in 1:length(vector_sp_2_tmp)){
              append(vector_sp_2,"_")
            }
            
            score <- spScore(vector_sp_1_tmp, vector_sp_2)
            score_table[z,x] <- score_table[z,x-1] - score
          }
          else if (x == 1 && z != 1 ){
            arrow_table[z,x] = "^"
            vector_sp_1 <- c()
            
            for (a in 1:length(vector_sp_1_tmp)){
              append(vector_sp_1,"_")
            }
            
            score <- spScore(vector_sp_1, vector_sp_2_tmp)
            score_table[z,x] <- score_table[z-1,x] - score
          }
          else {
            score <- spScore(vector_sp_1_tmp,vector_sp_2_tmp)
            
            vector_sp_1_gap <- c()
            vector_sp_2_gap <- c()
            
            for (a in 1:length(vector_sp_1_tmp)){
              vector_sp_1_gap <- append(vector_sp_1_gap,"_")
            }
            for (a in 1:length(vector_sp_2_tmp)){
              vector_sp_2_gap <- append(vector_sp_2_gap,"_")
            }
            
            gap_penalty_row <- spScore(vector_sp_1_tmp, vector_sp_2_gap)
            gap_penalty_col <- spScore(vector_sp_1_gap, vector_sp_2_tmp)
            
            if ( ( (score_table[z-1,x-1] + score  ) >= ( score_table[z,x-1] - gap_penalty_row ) ) && ((score_table[z-1,x-1] + score) >= (score_table[z-1,x] - gap_penalty_col ) ) ) {
              score_table[z,x] <- ( score_table[z-1,x-1] + score ) 
              arrow_table[z,x] <- "/"
              
            }
            else if ( ( (score_table[z-1,x] - gap_penalty_col ) >= ( score_table[z-1,x-1] + score ) ) && (score_table[z-1,x] - gap_penalty_col >= (score_table[z,x-1] - gap_penalty_row ) ) ) {
              score_table[z,x] <- ( score_table[z-1,x] - gap_penalty_col )
              arrow_table[z,x] <- "^"
            }
            else if (  ( (score_table[z,x-1] - gap_penalty_row ) >= ( score_table[z-1,x-1] + score) ) && (score_table[z,x-1] - gap_penalty_row >= (score_table[z-1,x] - gap_penalty_col ) ) ){
              score_table[z,x] <- ( score_table[z,x-1] - gap_penalty_row  )
              arrow_table[z,x] <- "<-"
            }
          }
          
        }
      }
      
      temp <- step[1]
      temp <- paste(temp , backTracking(arrow_table), sep ="/")
      
      alignments <- append(alignments, temp)
      #----------------------------------------------
    }
  }
  return(alignments)
}

printToTxt <- function(result){
  
  name <- strsplit(toString(result) , split = "/")[[1]][1]
  sequence <- strsplit(toString(result) , split = "/")[[1]][2]
  
  names <- strsplit(name , split = "-")[[1]]
  sequences <- strsplit(sequence , split = "-")[[1]]
  
  max_length <- 0
  
  for (i in 1 : length(names)){
    names[i] <- gsub(">" , "" , names[i])
    if (nchar(names[i] > max_length)){
      max_length <- nchar(names[i])
    }
  }
  
  lines <- c()
  
  for (i in 1 : length(names)){
    tmp <- ""
    chars_of_name <- strsplit(names[i], split ="")[[1]]
    
    for (j in 1:max_length){
      if (j <= length(chars_of_name)){
        tmp <- paste(tmp , chars_of_name[j], sep = "")
      }
      else {
        tmp <- paste(tmp , " ", sep = "")
      }
    }
    
    tmp <- paste(tmp , "    " , sep = "")
    tmp <- paste(tmp, sequences[i])
    lines <- append(lines, tmp)
  }
  
  fileConn <- file("MultipleAlignment.txt")
  writeLines(lines,fileConn)
  close(fileConn)
}

main <- function(){
  
  blosum62 <<- readBlosum62()
  gap_penalty <<- gapPenalty()
  print("PLEASE WAIT!")
  
  inputs<- readInput()
  
  seq_seq <- inputs[[1]]
  seq_names <- inputs[[2]]
  
  dst_mtrx <- allPairwiseAlignments(seq_seq, seq_names)
  guide_tree <- buildGuideTree(dst_mtrx)
  
  print("PLEASE KEEP WAITING")
  
  multipleAlignment <- progressiveAlignment(guide_tree , seq_names , seq_seq)
  
  printToTxt(multipleAlignment[length(multipleAlignment)])
  print("COMPLETED!")
  
}

main()