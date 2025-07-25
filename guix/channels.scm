(list
  (channel
    (name 'guix)
    (url "https://codeberg.org/guix/guix.git")
    (branch "master")
    (commit "216a37ba5005148bbb88c4f6b8e9dd5904d49074")
    (introduction
      (make-channel-introduction "9edb3f66fd807b096b48283debdcddccfea34bad"
        (openpgp-fingerprint "BBB0 2DDF 2CEA F6A8 0D1D  E643 A2A0 6DF2 A33A 54FA")
      )
    )
  )
  ; (channel ; somehow guix-past makes the build crash!, just excluding it works
  ;       (name 'guix-past)
  ;       (url "https://codeberg.org/guix-science/guix-past.git")
  ;       (branch "master")
  ;       (commit
  ;       "a6fc9859837de4a3d8cc824845f8bd28f68c530f")
  ;       (introduction
  ;       (make-channel-introduction
  ;       "0c119db2ea86a389769f4d2b9c6f5c41c027e336"
  ;       (openpgp-fingerprint
  ;       "3CE4 6455 8A84 FDC6 9DB4  0CFB 090B 1199 3D9A EBB5"))))
	(channel
    (name 'guix-cran)
    (url "https://github.com/guix-science/guix-cran.git")
    (branch "master")
    (commit "5c6128fc7c5a412c42eb165a4c641276d581b472")
  )

  (channel
    (name 'guix-science)
    (url "https://codeberg.org/guix-science/guix-science.git")
    (branch "master")
    (commit "31387efc4d90362e79251953fd1a6554dc262cda")
    (introduction
      (make-channel-introduction "b1fe5aaff3ab48e798a4cce02f0212bc91f423dc"
        (openpgp-fingerprint "CA4F 8CF4 37D7 478F DA05  5FD4 4213 7701 1A37 8446")
      )
    )
  )
    (channel
    (name 'guix-bioc)
    (url "https://github.com/guix-science/guix-bioc.git")
    (branch "master")
    (commit
        "700c5a6700fafd62b1fe4a789525b2053b1a236f")
    )
    
    ; (channel
    ;     (name 'guix-science-nonfree)
    ;     (url "https://codeberg.org/guix-science/guix-science-nonfree")
    ;     (branch "master")
    ;     (commit
    ;         "448d45309785272ba37510c6ef7b1dba9c4ddae4"))
)
			  
