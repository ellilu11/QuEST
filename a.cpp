nt i=0; i<iztrg.siztrge(); i++ ){
            if ((real(Ztrg[iztrg[i]]) <= real(center)) && (imag(Ztrg[iztrg[i]]) <= imag(center))){
                iztrgChild[0].push_back(iztrg[i]);
            } else if (real(Ztrg[iztrg[i]]) > real(center) && imag(Ztrg[iztrg[i]]) <= imag(center)){
                iztrgChild[1].push_back(iztrg[i]);
            } else if (real(Ztrg[iztrg[i]]) <= real(center) && imag(Ztrg[iztrg[i]]) > imag(center)){
                iztrgChild[2].push_back(iztrg[i]);
            } else if (real(Ztrg[iztrg[i]]) > real(center) && imag(Ztrg[iztrg[i]]) > imag(center)){
                iztrgChild[3].push_back(iztrg[i]);
            }
        }

