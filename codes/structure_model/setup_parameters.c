

void read_parameter_file(char *fname)
{
    printf("\nreading parameter file: %s!\n\n", fname);
    
    strcpy(tag[nt], "Mah_simu");
    add[nt] = &Mah_simu;
    id[nt++] = INT;
    
    strcpy(tag[nt], "Mah_dir");
    addr[nt] = Mah_dir;
    id[nt++] = STRING;
    
    strcpy(tag[nt], "Mah_file_basename");
    add[nt] = Mah_file_basename;
    id[nt++] = STRING;
    
    strcpy(tag[nt], "BaryonFrac");
    add[nt] = &BaryonFrac;
    id[nt++] = DOUBLE;
    
    strcpy(tag[nt], "Metal_gas_evolu");
    add[nt] = &Metal_gas_evolu;
    id[nt++] = INT;
    
    strcpy(tag[nt], "SolarMetallicity");
    add[nt] = &SolarMetallicity;
    id[nt++] = DOUBLE;
    
    strcpy(tag[nt], "MinimumMetallicityRelativeToSolar");
    add[nt] = &MinimumMetallicityRelativeToSolar;
    id[nt++] = DOUBLE;
    
    strcpy(tag[nt], "Do_preheating");
    add[nt] = &Do_preheating;
    id[nt++] = INT;
    
    strcpy(tag[nt], "Do_reinfall");
    add[nt] = &Do_reinfall;
    id[nt++] = INT;
    
    strcpy(tag[nt], "Redshift_end");
    add[nt] = &Redshift_end;
    id[nt++] = DOUBLE;
    


    
    
    if((fd = fopen(fname, "r")))
    {
        while(!feof(fd))
        {
            *buf = 0;
            fgets(buf, 200, fd);
            if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
                continue;
            
            if(buf1[0] == '%')
                continue;
            
            for(i = 0, j = -1; i < nt; i++)
                if(strcmp(buf1, tag[i]) == 0)
                {
                    j = i;
                    tag[i][0] = 0;
                    break;
                }
            
            if(j >= 0)
            {
                printf("%35s\t%10s\n", buf1, buf2);
                
                switch (id[j])
                {
                    case DOUBLE:
                        *((double *) addr[j]) = atof(buf2);
                        break;
                    case STRING:
                        strcpy(addr[j], buf2);
                        break;
                    case INT:
                        *((int *) addr[j]) = atoi(buf2);
                        break;
                }
            }
            else
            {
                printf("Error in file %s:   Tag '%s' not allowed or multiple defined.\n", fname, buf1);
                errorFlag = 1;
            }
        }
        fclose(fd);
        
        i = strlen(OutputDir);
        if(i > 0)
            if(OutputDir[i - 1] != '/')
                strcat(OutputDir, "/");
    }
    else
    {
        printf("Parameter file %s not found.\n", fname);
        errorFlag = 1;
    }

    for(i = 0; i < nt; i++)
    {
        if(*tag[i])
        {
            printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
            errorFlag = 1;
        }
    }
    
    if(errorFlag)
        exit(1);
}
	
