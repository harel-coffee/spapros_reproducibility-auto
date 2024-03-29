
Evaluation configurations

01 evaluation on Mad20              : Fig 1d,2e,S1,S2,S3,S5,S8a(lung)
02 evaluation on Mad20 (unified ct) : Fig S6b(Mad)
03 evaluation on Kra21 (unified ct) : Fig S6b(Kra)
04 evaluation on Mey22 (unified ct) : Fig S6b(Mey)
05 evaluation on Lit20              : Fig S8a(heart)
06 evaluation on Mey22, Airwaywall  : Fig S8b
07 evaluation on Asp19 sc           : Fig S8c(sc)
08 evaluation on Asp19 ST           : Fig S8c(ST)


Selections that were evaluated with each evaluation config

eval config | selection configs
------------|----------------------
     01     | 01, 02, 04, ext(lung)
     02     | 06a, 06b
     03     | 06a, 06b
     04     | 06a, 06b
     05     | 07, 08, ext(heart)
     06     | 09
     07     | 10b
     08     | 10b
     
reminder on selection run ids:
    01 basic feature selection methods lung (n150) : Fig 1d,S1,S2,S3,S8a(lung n150)
    02 spapros selection (n50 & 150)               : Fig 2c,2e(spapros),S4(without constraint),S5(spapros),S8a(lung)
    03 spapros selection with marker list          : Fig 2d
    04 basic feature selection methods lung (n50)  : Fig 2e(DE),S8a(lung n50)
    05 spapros selection with penalties            : Fig S4(with constraint)
    06a,06b cross data set selections              : Fig S6
    07 basic feature selection methods heart       : Fig S8a(heart)
    08 spapros selection heart                     : Fig S8a(heart)
    09 spapros selection airway wall               : Fig S8b
    10a,10b spapros selection dev heart            : Fig S8c
    ext = selections with external methods