citHeader("To cite package 'HiLDA' in publications use:")

year <- sub("-.*", "", meta$Date)
note <- sprintf("R package version %s", meta$Version)

citEntry(entry="Article",
         title = "HiLDA: a statistical approach to investigate differences in 
         mutational signatures",
         author = personList(as.person("Zhi", "Yang"),
                             as.person("Priyatama", "Pandey"),
                             as.person("Darryl", "Shibata"),
                             as.person("David V.", "Conti"),
                             as.person("Paul", "Marjoram"), 
                             as.person("Kimberly", "Siegmund")),
         year = 2019,
         journal = "bioRxiv",
         note = "(Currently under review for PeerJ)",
         url = "https://github.com/USCbiostats/HiLDA"
         )
         
