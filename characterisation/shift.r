# shift.r
# Shifts an array by a specified number of rows
# Used in computing MSDs etc.

shift = function(data,amount){
  
  # How many rows in input data
  rows <- nrow(data)
  
  # Split the array in two
  newbottom <- data[1:amount,]
  newtop <- data[(amount+1):rows,]
  
  # Recombine the two halves
  output <- rbind(newtop,newbottom)
  
  return(output)
  
  
}