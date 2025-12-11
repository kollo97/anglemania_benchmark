;; The following expression represents a list of channels.
(list(channel
        (name 'guix)
        (url "https://codeberg.org/guix/guix.git")
        (branch "master")
        (commit
          "364694773ef1c61d4c9ec73e24f73642786486ca")
        (introduction
          (make-channel-introduction
            "9edb3f66fd807b096b48283debdcddccfea34bad"
            (openpgp-fingerprint
              "BBB0 2DDF 2CEA F6A8 0D1D  E643 A2A0 6DF2 A33A 54FA"))))

      (channel
        (name 'guix-science)
        (url "https://codeberg.org/guix-science/guix-science.git")
        (branch "master")
        (commit
          "308d22d6ad9537be5684f5f6ba97833823fa7751"))
      (channel
        (name 'guix-cran)
        (url "https://github.com/guix-science/guix-cran.git")
        (branch "master")
        (commit
          "e6bedf1e7494001c5eef878a1c08ad8fdb11e6b7"))
      (channel
        (name 'guix-bioc)
        (url "https://github.com/guix-science/guix-bioc.git")
        (branch "master")
        (commit
          "c13272a395171ce8f76ad30f30da6d25f2c70c8b"))
)
