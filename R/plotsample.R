
plotsample <- function(a, n.layer=2){   
 
 if(n.layer==2) {

   par(mfrow=c(3,3))

   plot(a$expr, type="l", xlab="iteration", ylab="sample")
   title("expr")

   plot(a$mu, type="l", xlab="iteration",ylab="sample")
   title("mu")

   plot(a$mu, type="n", xlab="",ylab="",  xaxt="n", yaxt="n", bty="n")

   plot(a$gene, type="l", xlab="iteration",ylab="sample")
   title("gene")

   plot(a$cond, type="l", xlab="iteration",ylab="sample")
   title("cond")

   plot(a$inter, type="l", xlab="iteration",ylab="sample")
   title("inter")

   plot(a$var.e, type="l", xlab="iteration",ylab="sample")
   title("var.e")

   plot(a$var.b, type="l", xlab="iteration",ylab="sample")
   title("var.b")

 }


 if(n.layer==1) {

   par(mfrow=c(3,3))


   plot(a$mu, type="l", xlab="iteration",ylab="sample")
   title("mu")

   plot(a$gene, type="l", xlab="iteration",ylab="sample")
   title("gene")

   plot(a$cond, type="l", xlab="iteration",ylab="sample")
   title("cond")

   plot(a$inter, type="l", xlab="iteration",ylab="sample")
   title("inter")

   plot(a$var.t, type="l", xlab="iteration",ylab="sample")
   title("var.t")

 }

}
