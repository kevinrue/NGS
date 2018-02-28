
# Oscillate Y nudge above and below data points
# y: amount of Y nudge
# length: count of data points to nudge
nudge_balanced <- function(y, length){
  direction <- rep_len(c(-1,1), length.out=length)
  direction * y
}
