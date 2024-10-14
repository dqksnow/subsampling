(TeX-add-style-hook
 "references"
 (lambda ()
   (setq TeX-command-extra-options
         "-synctex=1")
   (LaTeX-add-bibitems
    "wang2021nonuniform"
    "fithian2014local"
    "wang2018optimal"
    "wang2019more"
    "wang2022maximum"
    "yao2019optimal"
    "yao2023model"))
 '(or :bibtex :latex))

