library(R.matlab)

    addresses=c("/VOT/","/Formantes/")
    
    for (ad in 1:2){
        load (file.path (paste("../Sujeito2",addresses[ad],"Identification/data/joao.dat",sep="")))

        beta <- fm$coefficients[2]

        filename <- paste("../Sujeito2",addresses[ad],"Identification/beta.mat",sep="")
        writeMat(filename,beta=beta)
    }

        
