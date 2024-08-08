#### Data Import ####


list.data <- lapply(name.file,
                    function(x) {
                      
                      tryCatch(
                        
                        {
                          
                          list(
                            # Stats Input
                            "Datasheet" = readxl::read_excel(x,
                                               sheet = 1),
                            
                            "Annotations" = tryCatch(
                              
                              {
                                
                                # Annotation Information
                                readxl::read_excel(x,
                                                   sheet = 2)
                                },
                              error = function(e) {
                                
                                print("Warning: Sheet 2 Not Detected")
                                }
                              )
                            )
                          },
                        
                        error = function(e) {
                          
                          print("Invalid Data Upload: Ensure Data is Included in Sheet 1 and Annotations in Sheet 2")
                          }
                        )
                      }
                    )

names(list.data) <- "Data Import"
