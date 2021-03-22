###Island & Group Name Harmonization function - by Michael Wohlwend, 12/20 ###
island.unifier<-function(island,islandgroup){
  isrep<-read.csv2("islands_rep_1.csv")
  
  #manually putting in phrases that frequently get read wrongly, depending on the users standard settings. It may not be pretty, but it is saver & easier to check
  colnames(isrep)[1]<-"Island_Phrase_Delet"#bug
  isrep$Island_Phrase_Delet[which(isrep$Island_Phrase_Delet=="\\\\s*\\\\([^\\\\)]+\\\\)")]<-"\\s*\\([^\\)]+\\)"
  isrep$Island_Asign_Submit[116]<-"San Crist�bal"
  isrep$Island_Asign_Return[116]<-"San Cristobal"
  isrep$Island_Asign_Submit[117]<-"San Cristóbal"
  isrep$Island_Asign_Return[117]<-"San Cristobal"
  isrep$Island_Asign_Submit[118]<-"Noreste de Santa F#"
  isrep$Island_Asign_Return[118]<-"Santa Fe"
  isrep$Island_Asign_Submit[119]<-"Pinz�n"
  isrep$Island_Asign_Return[119]<-"Pinzon"
  isrep$Island_Asign_Submit[120]<-"Pinzón"
  isrep$Island_Asign_Return[120]<-"Pinzon"
  isrep$Island_Asign_Submit[121]<-"R�bida"
  isrep$Island_Asign_Return[121]<-"Rabida"
  isrep$Island_Asign_Submit[122]<-"Santa Fé"
  isrep$Island_Asign_Return[122]<-"Santa Fe"
  isrep$Island_Asign_Submit[123]<-"Läna'i" 
  isrep$Island_Asign_Return[123]<-"Lana'i"
  isrep$Island_Asign_Submit[124]<-"Española" 
  isrep$Island_Asign_Return[124]<-"Espanola"
  islandgroup[which(islandgroup=="Juan Fernández Is.")]<-"Juan Fernandez"
  islandgroup[which(islandgroup=="Juan Fernández Is.")]<-"Juan Fernandez"
  islandgroup[which(islandgroup=="Galápagos")]<-"Galapagos"
  islandgroup[which(islandgroup=="San Crist�bal Islands")]<-"Galapagos"
  islandgroup[which(islandgroup=="Santa F# Islands")]<-"Galapagos"
  
  # Applying the list
  
  for(i in 1:nrow(isrep)){
    island<-gsub(isrep$Island_Phrase_Delet[i],"",island)
  } #detects pointless phrases (e.g. "Volacano" or ", Fijis") and non-included groups (e.g."Philipines") and replaces them with blanks
  for(i in which(island%in%isrep$Island_Asign_Submit)){
    island[i]<-unique(isrep$Island_Asign_Return[which(isrep$Island_Asign_Submit==island[i])])
  } # detects misspelled or otherwise differently named islands and harmonizes them. 
  for(i in which(island%in%isrep$Island_to_Group)){
    islandgroup[i]<-unique(isrep$Group_from_Island[which(isrep$Island_to_Group==island[i])])
    island[i]<-NA
  } #detects islands that are really groups and assigns them to the correct group
  for(i in which(islandgroup%in%isrep$Group_to_Group)){
    islandgroup[i]<-unique(isrep$Group_from_Group[which(isrep$Group_to_Group==islandgroup[i])])
  } #makes sure that every island is correctly placed. The step before is necessary because not all points have islands (some just groups)
  for(i in which(island%in%isrep$Island_final_Submit)){
    islandgroup[i]<-unique(isrep$Group_final_Asign[which(isrep$Island_final_Submit==island[i])])
  } #makes sure that every island is correctly placed. The step before is necessary because not all points have islands (some just groups)
  return_is<<-island
  return_isg<<-islandgroup
  isgrmat<-matrix(NA,length(island),2)
  isgrmat[,1]<-island
  isgrmat[,2]<-islandgroup
  colnames(isgrmat)<-c("island","islandgroup")
  return(isgrmat)
  }


