

AC_DEFUN([ABI_CC_OPTFLAGS],[
  dnl Init
  abi_cc_vendor_opt="none"
  abi_cc_version_opt="none"
  abi_cpu_spec_opt="none"

  dnl Look for optimizations
  AC_MSG_CHECKING([which cc optimizations to apply])

  dnl Case built from config/optim/cc_*.conf
  case "${abi_cc_vendor}" in
    ibm)
      abi_cc_vendor_opt="ibm"
      abi_cc_version_opt="default"
      case "${abi_cpu_spec}" in
        ibm_powerpc)
          abi_cpu_spec_opt="ibm_powerpc"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O4 -qarch=auto -qtune=auto -qstrict -qspill=2000"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2 -qarch=auto -qtune=auto -qstrict -qspill=2000"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O3 -qarch=auto -qtune=auto -qstrict -qspill=2000"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        ibm_powerpc64)
          abi_cpu_spec_opt="ibm_powerpc64"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O4 -qarch=auto -qtune=auto -qstrict -qspill=2000"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2 -qarch=auto -qtune=auto -qstrict -qspill=2000"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O3 -qarch=auto -qtune=auto -qstrict -qspill=2000"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    gnu)
      abi_cc_vendor_opt="gnu"
      abi_cc_version_opt="default"
      abi_cpu_spec_opt="default"
      case "${enable_optim}" in
        aggressive)
          enable_optim_opt="aggressive"
          CFLAGS_OPTIM="-O3 -mtune=native -march=native"
          ;;
        safe)
          enable_optim_opt="safe"
          CFLAGS_OPTIM="-O2"
          ;;
        standard)
          enable_optim_opt="standard"
          CFLAGS_OPTIM="-O2 -mtune=native -march=native"
          ;;
      esac   # [case: enable_optim, indent: 2, item: True]
      ;;
    open64)
      abi_cc_vendor_opt="open64"
      abi_cc_version_opt="default"
      case "${abi_cpu_spec}" in
        intel_pentium4)
          abi_cpu_spec_opt="intel_pentium4"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O3 -march=pentium4 -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2 -march=pentium4 -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O2 -march=pentium4 -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        amd_opteron)
          abi_cpu_spec_opt="amd_opteron"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O3 -march=opteron -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2 -march=opteron -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O2 -march=opteron -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    pathscale)
      abi_cc_vendor_opt="pathscale"
      abi_cc_version_opt="default"
      case "${abi_cpu_spec}" in
        intel_pentium4)
          abi_cpu_spec_opt="intel_pentium4"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O3 -march=pentium4 -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2 -march=pentium4 -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O2 -march=pentium4 -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        amd_opteron)
          abi_cpu_spec_opt="amd_opteron"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O3 -march=opteron -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2 -march=opteron -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O2 -march=opteron -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              CFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              CFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    intel)
      abi_cc_vendor_opt="intel"
      abi_cc_version_opt="default"
      abi_cpu_spec_opt="default"
      case "${enable_optim}" in
        aggressive)
          enable_optim_opt="aggressive"
          CFLAGS_OPTIM="-O3"
          ;;
        safe)
          enable_optim_opt="safe"
          CFLAGS_OPTIM="-O2"
          ;;
        standard)
          enable_optim_opt="standard"
          CFLAGS_OPTIM="-O2"
          ;;
      esac   # [case: enable_optim, indent: 2, item: True]
      ;;
  esac   # [case: abi_cc_vendor, indent: 0, item: True]

  dnl Display settings
  AC_MSG_RESULT([${abi_cc_vendor_opt}/${abi_cc_version_opt}/${abi_cpu_spec_opt}])

]) #ABI_CC_OPTFLAGS


AC_DEFUN([ABI_CXX_OPTFLAGS],[
  dnl Init
  abi_cxx_vendor_opt="none"
  abi_cxx_version_opt="none"
  abi_cpu_spec_opt="none"

  dnl Look for optimizations
  AC_MSG_CHECKING([which cxx optimizations to apply])

  dnl Case built from config/optim/cxx_*.conf
  case "${abi_cxx_vendor}" in
    ibm)
      abi_cxx_vendor_opt="ibm"
      abi_cxx_version_opt="default"
      case "${abi_cpu_spec}" in
        ibm_powerpc)
          abi_cpu_spec_opt="ibm_powerpc"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-O4 -qarch=auto -qtune=auto -qstrict -qspill=2000 -qessl"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-O2 -qarch=auto -qtune=auto -qstrict -qspill=2000 -qessl"
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-O3 -qarch=auto -qtune=auto -qstrict -qspill=2000 -qessl"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        ibm_powerpc64)
          abi_cpu_spec_opt="ibm_powerpc64"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-O4 -qarch=auto -qtune=auto -qstrict -qspill=2000 -qessl"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-O2 -qarch=auto -qtune=auto -qstrict -qspill=2000 -qessl"
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-O3 -qarch=auto -qtune=auto -qstrict -qspill=2000 -qessl"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    gnu)
      abi_cxx_vendor_opt="gnu"
      abi_cxx_version_opt="default"
      abi_cpu_spec_opt="default"
      case "${enable_optim}" in
        aggressive)
          enable_optim_opt="aggressive"
          CXXFLAGS_OPTIM="-O3 -mtune=native -march=native"
          ;;
        safe)
          enable_optim_opt="safe"
          CXXFLAGS_OPTIM="-O2"
          ;;
        standard)
          enable_optim_opt="standard"
          CXXFLAGS_OPTIM="-O2 -mtune=native -march=native"
          ;;
      esac   # [case: enable_optim, indent: 2, item: True]
      ;;
    pathscale)
      abi_cxx_vendor_opt="pathscale"
      abi_cxx_version_opt="default"
      case "${abi_cpu_spec}" in
        intel_pentium4)
          abi_cpu_spec_opt="intel_pentium4"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-O3 -march=pentium4 -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-O2 -march=pentium4 -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-O2 -march=pentium4 -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        amd_opteron)
          abi_cpu_spec_opt="amd_opteron"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-O3 -march=opteron -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-O2 -march=opteron -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-O2 -march=opteron -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              CXXFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              CXXFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              CXXFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    intel)
      abi_cxx_vendor_opt="intel"
      abi_cxx_version_opt="default"
      abi_cpu_spec_opt="default"
      case "${enable_optim}" in
        aggressive)
          enable_optim_opt="aggressive"
          CXXFLAGS_OPTIM="-O3"
          ;;
        safe)
          enable_optim_opt="safe"
          CXXFLAGS_OPTIM="-O2"
          ;;
        standard)
          enable_optim_opt="standard"
          CXXFLAGS_OPTIM="-O2"
          ;;
      esac   # [case: enable_optim, indent: 2, item: True]
      ;;
  esac   # [case: abi_cxx_vendor, indent: 0, item: True]

  dnl Display settings
  AC_MSG_RESULT([${abi_cxx_vendor_opt}/${abi_cxx_version_opt}/${abi_cpu_spec_opt}])

]) #ABI_CXX_OPTFLAGS


AC_DEFUN([ABI_FC_OPTFLAGS],[
  dnl Init
  abi_fc_vendor_opt="none"
  abi_fc_version_opt="none"
  abi_cpu_spec_opt="none"

  dnl Look for optimizations
  AC_MSG_CHECKING([which fc optimizations to apply])

  dnl Case built from config/optim/fc_*.conf
  case "${abi_fc_vendor}" in
    ibm)
      abi_fc_vendor_opt="ibm"
      abi_fc_version_opt="default"
      case "${abi_cpu_spec}" in
        ibm_powerpc)
          abi_cpu_spec_opt="ibm_powerpc"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O4 -qmaxmem=65536 -qspill=2000 -qarch=auto -qtune=auto -qcache=auto -qstrict -qsuppress=1520-022:1520-031"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -qmaxmem=65536 -qspill=2000 -qarch=auto -qtune=auto -qcache=auto -qstrict -qsuppress=1520-022:1520-031"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O3 -qmaxmem=65536 -qspill=2000 -qarch=auto -qtune=auto -qcache=auto -qstrict -qsuppress=1520-022:1520-031"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        ibm_powerpc64)
          abi_cpu_spec_opt="ibm_powerpc64"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O4 -qmaxmem=65536 -qspill=2000 -qarch=auto -qtune=auto -qcache=auto -qstrict -qsuppress=1520-022:1520-031"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -qmaxmem=65536 -qspill=2000 -qarch=auto -qtune=auto -qcache=auto -qstrict -qsuppress=1520-022:1520-031"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O3 -qmaxmem=65536 -qspill=2000 -qarch=auto -qtune=auto -qcache=auto -qstrict -qsuppress=1520-022:1520-031"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    gnu)
      abi_fc_vendor_opt="gnu"
      case "${abi_fc_version}" in
        4.7)
          abi_fc_version_opt="4.7"
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -mtune=native -march=native -faggressive-function-elimination -fstack-arrays"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -mtune=native -march=native"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_fc_version_opt="default"
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -mtune=native -march=native -funroll-loops -faggressive-function-elimination"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -mtune=native -march=native"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_fc_version, indent: 2, item: True]
      ;;
    open64)
      abi_fc_vendor_opt="open64"
      abi_fc_version_opt="default"
      case "${abi_cpu_spec}" in
        intel_pentium4)
          abi_cpu_spec_opt="intel_pentium4"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=pentium4 -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=pentium4 -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=pentium4 -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        amd_opteron)
          abi_cpu_spec_opt="amd_opteron"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=opteron -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=opteron -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=opteron -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -OPT:Olimit=0 -g -ggdbm"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    pgi)
      abi_fc_vendor_opt="pgi"
      abi_fc_version_opt="default"
      abi_cpu_spec_opt="default"
      case "${enable_optim}" in
        aggressive)
          enable_optim_opt="aggressive"
          FCFLAGS_OPTIM="-O3 nehalem-64 -m64"
          ;;
        safe)
          enable_optim_opt="safe"
          FCFLAGS_OPTIM="-O2 -tp nehalem-64 -m64"
          ;;
        standard)
          enable_optim_opt="standard"
          FCFLAGS_OPTIM="-O2 -tp nehalem-64 -m64"
          ;;
      esac   # [case: enable_optim, indent: 2, item: True]
      ;;
    g95)
      abi_fc_vendor_opt="g95"
      abi_fc_version_opt="default"
      case "${abi_cpu_spec}" in
        intel_pentium4)
          abi_cpu_spec_opt="intel_pentium4"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=pentium4 -mmmx -msse -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=pentium4 -mmmx -msse -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=pentium4 -mmmx -msse -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_coreduo)
          abi_cpu_spec_opt="intel_coreduo"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=prescott -mmmx -msse -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=prescott -mmmx -msse -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=prescott -mmmx -msse -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_xeon)
          abi_cpu_spec_opt="intel_xeon"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=nocona -mmmx -msse -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=nocona -mmmx -msse"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O3 -march=nocona -mmmx -msse -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        amd_opteron)
          abi_cpu_spec_opt="amd_opteron"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=opteron"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=opteron"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=opteron"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_itanium1)
          abi_cpu_spec_opt="intel_itanium1"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O3"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        dec_alphaev67)
          abi_cpu_spec_opt="dec_alphaev67"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O3"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        ibm_powerpc)
          abi_cpu_spec_opt="ibm_powerpc"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O4 -mpowerpc"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -mpowerpc"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -mpowerpc"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        dec_alphaev56)
          abi_cpu_spec_opt="dec_alphaev56"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O3"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_itanium2)
          abi_cpu_spec_opt="intel_itanium2"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O3"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        amd_athlon)
          abi_cpu_spec_opt="amd_athlon"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=athlon"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=athlon"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=athlon"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_core2)
          abi_cpu_spec_opt="intel_core2"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=prescott -mmmx -msse -msse2 -msse3"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=prescott -mmmx -msse -msse2 -msse3"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=prescott -mmmx -msse -msse2 -msse3"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        amd_athlon64)
          abi_cpu_spec_opt="amd_athlon64"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=athlon64"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=athlon64"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=athlon64"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_pentium3)
          abi_cpu_spec_opt="intel_pentium3"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=pentium3 -mmmx -msse"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=pentium3 -mmmx -msse"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=pentium3 -mmmx -msse"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        ibm_powerpc64)
          abi_cpu_spec_opt="ibm_powerpc64"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O4 -mpowerpc64"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -mpowerpc64"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -mpowerpc64"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    fujitsu)
      abi_fc_vendor_opt="fujitsu"
      abi_fc_version_opt="default"
      abi_cpu_spec_opt="default"
      case "${enable_optim}" in
        aggressive)
          enable_optim_opt="aggressive"
          FCFLAGS_OPTIM="-Of -X9 -Ps -Wv,-md"
          ;;
        safe)
          enable_optim_opt="safe"
          FCFLAGS_OPTIM="-Of -X9 -Ps -Wv,-md"
          ;;
        standard)
          enable_optim_opt="standard"
          FCFLAGS_OPTIM="-Of -X9 -Ps -Wv,-md"
          ;;
      esac   # [case: enable_optim, indent: 2, item: True]
      ;;
    pathscale)
      abi_fc_vendor_opt="pathscale"
      abi_fc_version_opt="default"
      case "${abi_cpu_spec}" in
        intel_pentium4)
          abi_cpu_spec_opt="intel_pentium4"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=pentium4 -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=pentium4 -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=pentium4 -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        amd_opteron)
          abi_cpu_spec_opt="amd_opteron"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -march=opteron -msse2"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -march=opteron -msse2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -march=opteron -msse2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
    nag)
      abi_fc_vendor_opt="nag"
      abi_fc_version_opt="default"
      abi_cpu_spec_opt="default"
      case "${enable_optim}" in
        aggressive)
          enable_optim_opt="aggressive"
          FCFLAGS_OPTIM="-O4"
          ;;
        safe)
          enable_optim_opt="safe"
          FCFLAGS_OPTIM="-O2"
          ;;
        standard)
          enable_optim_opt="standard"
          FCFLAGS_OPTIM="-O3"
          ;;
      esac   # [case: enable_optim, indent: 2, item: True]
      ;;
    intel)
      abi_fc_vendor_opt="intel"
      abi_fc_version_opt="default"
      case "${abi_cpu_spec}" in
        intel_itanium1)
          abi_cpu_spec_opt="intel_itanium1"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -fp-model fast=1 -fp-relaxed -ip"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -fp-model precise -fp-speculation=safe"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        intel_itanium2)
          abi_cpu_spec_opt="intel_itanium2"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -fp-model fast=1 -fp-relaxed -ip"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -fp-model precise -fp-speculation=safe"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        amd_opteron)
          abi_cpu_spec_opt="amd_opteron"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O1"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
        *)
          abi_cpu_spec_opt="default"
          case "${enable_optim}" in
            aggressive)
              enable_optim_opt="aggressive"
              FCFLAGS_OPTIM="-O3 -xHOST"
              ;;
            safe)
              enable_optim_opt="safe"
              FCFLAGS_OPTIM="-O2 -xHost -fltconsistency -fp-model precise -fp-speculation=safe -prec-div -prec-sqrt"
              ;;
            standard)
              enable_optim_opt="standard"
              FCFLAGS_OPTIM="-O2 -xHost"
              ;;
          esac   # [case: enable_optim, indent: 4, item: True]
          ;;
      esac   # [case: abi_cpu_spec, indent: 2, item: True]
      ;;
  esac   # [case: abi_fc_vendor, indent: 0, item: True]

  dnl Display settings
  AC_MSG_RESULT([${abi_fc_vendor_opt}/${abi_fc_version_opt}/${abi_cpu_spec_opt}])

]) #ABI_FC_OPTFLAGS
