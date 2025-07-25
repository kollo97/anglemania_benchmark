;; What follows is a "manifest" equivalent to the command line you gave.
;; You can store it in a file that you may then pass to any 'guix' command
;; that accepts a '--manifest' (or '-m') option.
(use-modules (gnu packages)
             (guix packages)
             ((guix licenses) #:prefix license:)
             (guix git-download)
             (guix build-system r)
             (guix download)
             (guix build-system python)
             (guix build-system pyproject))


(define-public r-anndatar
  (let ((commit "8a53b29102038fd9692cba02aa01b276f5bd58f0")  ;; Replace with the actual commit hash for the desired version
        (revision "1"))
    (package
      (name "r-anndatar")
      (version (git-version "0.99.0" revision commit))  ;; Replace with the correct version
      (source
       (origin
         (method git-fetch)
         (uri (git-reference
               (url "https://github.com/scverse/anndataR")  ;; Repository URL
               (commit commit)))
         (file-name (git-file-name name version))
         (sha256
          (base32 "0fp29mki2qk0pd5cgg6b1jn6ixwr3s17ny0fyi8z5qb75gxa9pi5"))))  ;; Replace with the actual SHA256 checksum
      (build-system r-build-system)
      (native-inputs (map specification->package
                          (list "r-knitr"              ;; from Suggests
                                "r-anndata"
                                "r-biocstyle"
                                "r-reticulate"          ;; >= 1.36.1
                                "r-hdf5r"
                                "r-matrix"               ;; >= 1.3.11
                                "r-rmarkdown"
                                "r-s4vectors"
                                "r-seuratobject"
                                "r-singlecellexperiment"
                                "r-summarizedexperiment"
                                "r-testthat"            ;; >= 3.0.0
                                "r-withr"
                                "r-tidyverse")))
      (home-page "https://github.com/theislab/anndataR")
      (synopsis "R interface for AnnData objects")
      (description
       "anndataR provides an R interface to work with AnnData objects, allowing integration with single-cell analysis pipelines in R.")
      (license license:gpl3))))

(define-public python-plottable
  (package
    (name "python-plottable")
    (version "0.1.5")
    (source
     (origin
       (method url-fetch)
       (uri (pypi-uri "plottable" version))
       (sha256
        (base32 "0784hkl5vgii0b623yfhjd1iwjx10g0hahppbgf2j8f864m7cp93"))))
    (build-system pyproject-build-system)
    (propagated-inputs (map specification->package
                              (list "python-setuptools"
                                    "python-wheel"
                                    "python-matplotlib"
                                    "python-numpy"
                                    "python-pandas"
                                    "python-pillow")))
    (native-inputs (map specification->package
                        (list "python-black"
                              "python-pytest")))
    (home-page "https://github.com/znstrider/plottable")
    (synopsis "Beautifully customized tables with matplotlib")
    (description "Beautifully customized tables with matplotlib.")
    (license license:expat)))

(define-public python-svgwrite
  (package
    (name "python-svgwrite")
    (version "1.4.3")
    (source
     (origin
       (method url-fetch)
       (uri (pypi-uri "svgwrite" version ".zip"))
       (sha256
        (base32 "1hr4wn68ph0yi5y776wdcbrdm0y6gy9vqxkzk9hsc0ik8kadzyx8"))))
    (build-system pyproject-build-system)
    (arguments
     `(#:tests? #f)) ; Disable the 'check' phase
    (propagated-inputs (map specification->package
                        (list "python-setuptools"
                              "python-wheel"
                              "zip"
                              "unzip")))
    (home-page "http://github.com/mozman/svgwrite.git")
    (synopsis "A Python library to create SVG drawings.")
    (description
     "This package provides a Python library to create SVG drawings.")
    (license license:expat)))

(define-public python-tree
  (package
    (name "python-tree")
    (version "0.2.4")
    (source
     (origin
       (method url-fetch)
       (uri (pypi-uri "Tree" version))
       (sha256
        (base32 "10vqdxj2gpns1iciym86fx78c46i4dd94y6sa7snkpahpz4qwkgq"))))
    (build-system pyproject-build-system)
    (propagated-inputs (cons* python-svgwrite
                        (map specification->package
                              (list "python-click"
                                    "python-pillow"
                                    "python-wheel"
                                    "python-setuptools"))))
    (home-page "https://github.com/PixelwarStudio/PyTree")
    (synopsis "A package for creating and drawing trees")
    (description
     "This package provides a package for creating and drawing trees.")
    (license license:expat)))

(define-public python-scib-metrics
  (package
    (name "python-scib-metrics")
    (version "0.5.1")
    (source
     (origin
       (method url-fetch)
       (uri (pypi-uri "scib_metrics" version))
       (sha256
        (base32 "0rp2l4rd44ghrvv7iyshg4fqh5skwmiz1yllp4119hgimi8h5lbl"))))
    (build-system pyproject-build-system)
    (arguments
     `(#:tests? #f  ; Disable the test phase
       #:phases
       (modify-phases %standard-phases
         ;; Add a phase to patch the license field in pyproject.toml
         (add-before 'build 'fix-license
           (lambda _
             (substitute* "pyproject.toml"
               (("\"license\"\\s*=\\s*\"\"") "\"license\" = \"MIT\""))
             #t))
         ;; Prevent installation of optional dependencies like 'scib'
         (add-before 'build 'set-hatch-env
           (lambda _
             (setenv "HATCH_EXTRAS_REQUIRED" "")
             #t)))))
    (propagated-inputs (cons* 
                              python-tree
                              python-plottable
                         (map specification->package
                            (list "python-anndata"
                                  "python-chex"
                                  "python-dm-tree"
                                  "python-igraph"
                                  "python-jax"
                                  "python-jaxlib"
                                  "python-matplotlib"
                                  "python-numpy"
                                  "python-pandas"
                                  "python-pynndescent"
                                  "python-rich"
                                  "python-scanpy"
                                  "python-scikit-learn"
                                  "python-scipy"
                                  "python-tqdm"
                                  "python-treelib"
                                  "python-umap-learn"
                                  "python-hatchling"))))
    (native-inputs (map specification->package
                        (list "python-wheel"
                              "python-setuptools"
                              "python-black"
                              "python-bump2version"
                              "python-coverage"
                              "python-flake8"
                              "python-joblib"
                              "python-numba"
                              "pre-commit"
                              "python-pytest"
                              "python-twine")))
    (home-page "https://github.com/yoseflab/scib-metrics")
    (synopsis "Accelerated and Python-only scIB metrics")
    (description "Accelerated and Python-only scIB metrics.")
    (license license:expat)))

(define-public python-scib
  (package
    (name "python-scib")
    (version "1.1.7")
    (source
     (origin
       (method url-fetch)
       ;; Assuming VERSION.txt is bundled in PyPI sdist
       (uri (pypi-uri "scib" version))
       (sha256
        (base32 "00h2k6cfnp3g6m3ap3xhfi2sk2j1x5rimfkv65f2dpwsp3bgxm9v"))))
    (build-system pyproject-build-system)
    (arguments
     `(#:tests? #f ; Disable tests by default, since they require extras
       #:phases
       (modify-phases %standard-phases
         (delete 'sanity-check)))) 
    (propagated-inputs
     (map specification->package
          (list "python-numpy"
                "python-pandas"
                "python-seaborn"
                "python-matplotlib"
                "python-numba"
                "python-scanpy"
                "python-anndata"
                "python-h5py"
                "python-scipy"
                "python-scikit-learn"
                "python-scikit-misc"
                "python-leidenalg"
                "python-umap-learn"
                "python-pydot"
                "python-igraph"
                "python-llvmlite"
                "python-deprecated")))
    (native-inputs
     (map specification->package
          (list "python-setuptools"
                "python-wheel")))
    (home-page "https://github.com/theislab/scib")
    (synopsis "Evaluating single-cell data integration methods")
    (description "scIB provides benchmarking metrics for evaluating single-cell data integration methods.")
    (license license:expat)))


(define-public r-anglemania
  (let ((commit "8052de59deae277d8941493bd72e026f04c5d854")  ;; Replace with the actual commit hash for the desired version
        (revision "1"))
    (package
      (name "r-anglemania")
      (version (git-version "0.99.0" revision commit))  ;; Replace with the correct version
      (source
       (origin
         (method git-fetch)
         (uri (git-reference
               (url "https://github.com/BIMSBbioinfo/anglemania/")  ;; Repository URL
               (commit commit)))
         (file-name (git-file-name name version))
         (sha256
          (base32 "1cg583wfacd259qjqizzs7l49f39z71kl5g9y553a708rsw8j61v"))))  ;; Replace with the actual SHA256 checksum
      (build-system r-build-system)
      ;; Disable tests here:
      (arguments
      '(#:tests? #f))
      (native-inputs (map specification->package
                          (list "r-matrix"
                                "r-checkmate"
                                "r-seurat"
                                "r-seuratobject"
                                "r-s4vectors"
                                "r-singlecellexperiment"
                                "r-summarizedexperiment"
                                "r-bigstatsr"
                                "r-dplyr"
                                "r-magrittr"
                                "r-tidyr"
                                "r-pbapply"
                                "r-testthat"
                                "r-withr"
                                "r-tidyverse"
                                "r-digest"
                                "r-matrixstats")))
      (home-page "https://github.com/BIMSBbioinfo/anglemania")
      (synopsis "New approach to the integration of scRNA-seq")
      (description
       "anglemania is a new approach to the integration of scRNA-seq and, potentially, others sc-omics from similar biological entities. The novelty, as well as the cornerstone, of the proposed approach, is to use the conservation of angles between gene pairs across an assembly of datasets to be integrated.")
      (license license:gpl3))))



(define python-stuff
  (list "python"
        "python-scanpy"
        "python-loompy"
        "python-pynvim"
        "python-anndata"
        "python-pandas"
        "python-dm-tree"
        "python-treelib"
        "python-chex"
        "python-scanorama"
        ; "python-pytorch-lightning"
        "python-scvi-tools"
        "python-harmonypy"
        "python-ipykernel"
        ))

(define r-stuff
  (list "r-minimal"
        "r-data-table" 
        "r-optparse"
        "r-purrr"
        "r-dplyr"
        "r-ggplot2"
        "r-tidyr"
        "r-ggbreak"
        "r-jsonlite"
        "r-r-utils"
        "r-unix"
        "r-rcpp"
        "r-stringr"
        "r-pbapply"
        "r-glmgampoi"
        "r-rmarkdown"
        "r-knitr"
        "r-languageserver"
        "r-devtools"
        "r-bigstatsr"
        "r-seurat"
        "r-rappdirs"
        "r-httpgd"
        "r-reticulate"
        "r-countsplit"
        "r-singlecellexperiment"
        "r-anndata"
        "r-tidyverse"
        "r-janitor"
        "r-hdf5r"
        "r-rcolorbrewer"
        "r-viridis"
        "r-uwot"
        "r-complexheatmap"
        "r-interactivecomplexheatmap"
        "r-magick"
        "r-kbet"
        "r-checkmate"
        "r-splatter"
        "r-jsonlite"
        "r-domc"
        "r-digest"
        "r-corrplot"
        ))

(define other-stuff
  (list "imagemagick"
        "zip"
        "snakemake"
        "guix-jupyter"
  ))


(packages->manifest
 (cons* r-anndatar
        python-plottable
        python-svgwrite
        python-tree
        python-scib-metrics
        python-scib
        ; r-anglemania
        (map specification->package
             (append python-stuff
                     r-stuff
                     other-stuff))))