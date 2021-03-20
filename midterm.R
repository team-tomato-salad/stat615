# midterm proj
library(plotly)

### load exmaple
# rm(list = ls())

images = load("~/Desktop/615/midterm/images_examples.RData")

# ??? image size will change
### Function starts here
midterm_proj = function(noisy_img){
  
n = dim(noisy_img)[1]
inside = n-1
### denoise part 
center_one = min(noisy_img)
center_two = max(noisy_img)

# four courners
noisy_img[1,1] =(((noisy_img[1,2]+noisy_img[2,1])/2)+noisy_img[1,1])/2
noisy_img[n,n] = (((noisy_img[inside,n]+noisy_img[n,inside])/2)+noisy_img[n,n])/2
noisy_img[1,n] =  (((noisy_img[1,inside]+noisy_img[2,n])/2)+noisy_img[1,n])/2
noisy_img[n,1] =  (((noisy_img[inside,1]+noisy_img[n,2])/2)+noisy_img[n,1])/2
# edge cases
for (e in 2:inside){
  forw = e-1
  nex = e+1
  noisy_img[1,e] =  (((noisy_img[1,forw]+noisy_img[1,nex]+noisy_img[2,e])/3)+noisy_img[1,e])/2
  noisy_img[n,e] =(((noisy_img[n,forw]+noisy_img[n,nex]+noisy_img[inside,e])/3)+noisy_img[1,e])/2
  noisy_img[e,1] =(((noisy_img[forw,1]+noisy_img[nex,1]+noisy_img[e,2])/3)+noisy_img[e,1])/2
  noisy_img[e,n] = (((noisy_img[forw,n]+noisy_img[nex,n]+noisy_img[e,inside])/3)+noisy_img[e,n])/2
}

# every pixel in between
for (i in 2:inside){
  # i is row index
  for (j in 2:inside){
    neighbours = (noisy_img[i+1,j]+noisy_img[i,j+1]+noisy_img[i-1,j]+noisy_img[i,j-1])/4
    noisy_img[i,j] = 0.5*noisy_img[i,j]+0.5*neighbours
      #(noisy_img[i,j]+noisy_img[i+1,j]+noisy_img[i,j+1]+noisy_img[i-1,j]+noisy_img[i,j-1]+noisy_img[i+1,j+1]+noisy_img[i-1,j-1])/5
  }
}

 plot_ly(z = noisy_img,type = "heatmap",colors = "Blues") %>%   
  layout(title ="Patch fixed image")

 
### k_means func

kmean = function(m1,m2,martix){
  diff = 2
  lap = 0
  #for (i in  1:5){
  while(diff>0.0005 & lap < 2){
    dist_to_one = abs(martix-m1)
    dist_to_two = abs(martix-m2)
    
    class_one = (dist_to_one-dist_to_two)>=0
    class_two = (dist_to_one-dist_to_two)<0
    
    one = martix[ which(class_one == T)]
    two = martix[ which(class_two== T)]
    past1 = m1
    past2 = m2
    m1 = mean(one)
    m2 = mean(two)
    sigma = round(sd(one))
    
    matrix = m1*class_one+m2*class_two
    diff = abs(past1+past2-m1-m2)
    lap = lap +1  
  }
  
  m1 = round(m1)
  m2 = round(m2)
  matrix = m1*class_one+m2*class_two
  # print(m1)
  # print(m2)
  # print(lap)
  return (list(mu1 = m1,mu2 = m2,ann = matrix))
}

result = kmean(center_one,center_two,noisy_img)

denoi = result$ann

### final rescue
for (i in 2:inside){
  # i is row index
  for (j in 2:inside){
    
    if(denoi[i+1,j]!=denoi[i,j] && denoi[i,j+1]!= denoi[i,j] && denoi[i-1,j]!= denoi[i,j] &&denoi[i,j-1]!=denoi[i,j]){
      denoi[i,j] = denoi[i+1,j]
    }
  }
}
  return (list(final_matrix = denoi,m1 = result$mu1,m2=result$mu2))
} # end of midterm_proj function 


#### run package here
denoi_mp = midterm_proj(noisy_img)
print(denoi_mp$m1)
print(denoi_mp$m2)

#####plotting result and calculate error 
plot_ly(z = denoi_mp$final_matrix,type = "heatmap",colors = "Blues") %>%   
  layout(title ="Kmeans after denoise")

sum(abs( denoi_mp$final_matrix-true_img))/(500^2)





### self_made example circle
m = matrix(1:(400^2),400,400)

g= expand.grid(1:nrow(m), 1:nrow(m))
g$d2 = sqrt ((g$Var1-280)^2 + (g$Var2-250)^2)
g$inside = g$d2<=50

m[which (g$inside == T)]=2
m[which (g$inside == F)]=10

tru = m
plot_ly(z = tru,type = "heatmap",colors = "Blues") %>%   
  layout(title ="test_tru")

noi = rnorm(400^2,sd = 2)
dim(noi) = c(400,400)

noi = tru+noi
plot_ly(z = noi,type = "heatmap",colors = "Blues") %>%   
  layout(title ="test_bur")
noisy_img = noi

### self_made example band

m= rep(40,500^2)
dim(m) = c(500,500)
m[12:130,] = 45
tru = m
plot_ly(z = tru,type = "heatmap",colors = "Blues") %>%   
  layout(title ="test_tru")
noi = rnorm(500^2,sd = 2)
dim(noi) = c(500,500)

noi = tru+noi
plot_ly(z = noi,type = "heatmap",colors = "Blues") %>%   
  layout(title ="test_bur")
noisy_img = noi

### self_made example square

m= rep(45,100^2)
dim(m) = c(100,100)
m[(15:30),(25:70)] = 20
m[(30:60),(25:45)] = 20
m[(60:70),(25:70)] = 20
tru = m
plot_ly(z = tru,type = "heatmap",colors = "Blues") %>%layout(title ="test_tru")
noi = rnorm(100^2,sd = 5)
dim(noi) = c(100,100)

noi = tru+noi
plot_ly(z = noi,type = "heatmap",colors = "Blues") %>%   
  layout(title ="test_bur")

noisy_img = noi

### self_made example square
another  = true_img
another[which (true_img == 3)]=0
another[which (true_img == 0)]=3

n = matrix(rnorm(500^2),500,500)
noisy_img = another+n
