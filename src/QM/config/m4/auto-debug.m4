

AC_DEFUN([ABI_CC_DBGFLAGS],[
  dnl Init
  abi_cc_vendor_dbg="none"
  abi_cc_version_dbg="none"
  abi_cpu_spec_dbg="none"

  dnl Look for debug flags
  AC_MSG_CHECKING([which cc debug flags to apply])

  dnl Case built from config/debug/cc_*.conf
  if test "${abi_cc_vendor}" = "gnu"; then
    abi_cc_vendor_dbg="gnu"
    case "${abi_cc_version}" in
      4.9)
        abi_cc_version_dbg="4.9"
        abi_cpu_spec_dbg="default"
        case "${enable_debug}" in
          naughty)
            enable_debug_dbg="naughty"
            CFLAGS_DEBUG="${CFLAGS_DEBUG} -g3 -ggdb -fdiagnostics-color=auto -Wall -Wextra -fbounds-check"
            ;;
          verbose)
            enable_debug_dbg="verbose"
            
            ;;
          paranoid)
            enable_debug_dbg="paranoid"
            CFLAGS_DEBUG="${CFLAGS_DEBUG} -g3 -ggdb -fdiagnostics-color=auto -Wall -Wextra"
            ;;
          enhanced)
            enable_debug_dbg="enhanced"
            CFLAGS_DEBUG="${CFLAGS_DEBUG} -g3 -ggdb -fdiagnostics-color=auto"
            ;;
          basic)
            enable_debug_dbg="basic"
            
            ;;
        esac   # [case: enable_debug, indent: 3, item: True]
        ;;
      *)
        abi_cc_version_dbg="default"
        abi_cpu_spec_dbg="default"
        case "${enable_debug}" in
          naughty)
            enable_debug_dbg="naughty"
            CFLAGS_DEBUG="${CFLAGS_DEBUG} -g3 -ggdb -Wall -Wextra -fbounds-check"
            ;;
          verbose)
            enable_debug_dbg="verbose"
            
            ;;
          paranoid)
            enable_debug_dbg="paranoid"
            CFLAGS_DEBUG="${CFLAGS_DEBUG} -g3 -ggdb -Wall -Wextra"
            ;;
          enhanced)
            enable_debug_dbg="enhanced"
            CFLAGS_DEBUG="${CFLAGS_DEBUG} -g3 -ggdb"
            ;;
          basic)
            enable_debug_dbg="basic"
            
            ;;
        esac   # [case: enable_debug, indent: 3, item: True]
        ;;
    esac   # [case: abi_cc_version, indent: 1, item: False]
  fi

  dnl Display settings
  AC_MSG_RESULT([${abi_cc_vendor_dbg}/${abi_cc_version_dbg}/${abi_cpu_spec_dbg}])

]) #ABI_CC_DBGFLAGS


AC_DEFUN([ABI_CXX_DBGFLAGS],[
  dnl Init
  abi_cxx_vendor_dbg="none"
  abi_cxx_version_dbg="none"
  abi_cpu_spec_dbg="none"

  dnl Look for debug flags
  AC_MSG_CHECKING([which cxx debug flags to apply])

  dnl WARNING: no config files were found for language

  dnl Display settings
  AC_MSG_RESULT([${abi_cxx_vendor_dbg}/${abi_cxx_version_dbg}/${abi_cpu_spec_dbg}])

]) #ABI_CXX_DBGFLAGS


AC_DEFUN([ABI_FC_DBGFLAGS],[
  dnl Init
  abi_fc_vendor_dbg="none"
  abi_fc_version_dbg="none"
  abi_cpu_spec_dbg="none"

  dnl Look for debug flags
  AC_MSG_CHECKING([which fc debug flags to apply])

  dnl Case built from config/debug/fc_*.conf
  case "${abi_fc_vendor}" in
    ibm)
      abi_fc_vendor_dbg="ibm"
      abi_fc_version_dbg="default"
      abi_cpu_spec_dbg="default"
      case "${enable_debug}" in
        naughty)
          enable_debug_dbg="naughty"
          FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -qnooptimize -qextcheck -qflag=i:i -qfloat=nans -qinitauto=7FBFFFFF -qflttrap=overflow:underflow:zerodivide:invalid:enable -qsigtrap -C -qcheck"
          ;;
        verbose)
          enable_debug_dbg="verbose"
          
          ;;
        paranoid)
          enable_debug_dbg="paranoid"
          FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -qnooptimize -qextcheck -qflag=i:i -qfloat=nans -qinitauto=7FBFFFFF -qflttrap=overflow:underflow:zerodivide:invalid:enable -qsigtrap"
          ;;
        enhanced)
          enable_debug_dbg="enhanced"
          FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -qnooptimize -qextcheck -qflag=i:i -qfloat=nans -qinitauto=7FBFFFFF"
          ;;
        basic)
          enable_debug_dbg="basic"
          
          ;;
      esac   # [case: enable_debug, indent: 2, item: True]
      ;;
    gnu)
      abi_fc_vendor_dbg="gnu"
      case "${abi_fc_version}" in
        4.9)
          abi_fc_version_dbg="4.9"
          abi_cpu_spec_dbg="default"
          case "${enable_debug}" in
            naughty)
              enable_debug_dbg="naughty"
              FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -fbacktrace -finit-real=nan -Wimplicit-interface -Wtabs -Wsurprising -fdiagnostics-color=auto -Wall -Wextra -Wfunction-elimination -fbounds-check"
              ;;
            verbose)
              enable_debug_dbg="verbose"
              
              ;;
            paranoid)
              enable_debug_dbg="paranoid"
              FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -fbacktrace -finit-real=nan -Wimplicit-interface -Wtabs -Wsurprising -fdiagnostics-color=auto -Wall -Wextra -Wfunction-elimination"
              ;;
            enhanced)
              enable_debug_dbg="enhanced"
              FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -fbacktrace -finit-real=nan -Wimplicit-interface -Wtabs -Wsurprising -fdiagnostics-color=auto"
              ;;
            basic)
              enable_debug_dbg="basic"
              
              ;;
          esac   # [case: enable_debug, indent: 4, item: True]
          ;;
        *)
          abi_fc_version_dbg="default"
          abi_cpu_spec_dbg="default"
          case "${enable_debug}" in
            naughty)
              enable_debug_dbg="naughty"
              FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -fbacktrace -finit-real=nan -Wimplicit-interface -Wtabs -Wall -Wextra -fbounds-check"
              ;;
            verbose)
              enable_debug_dbg="verbose"
              
              ;;
            paranoid)
              enable_debug_dbg="paranoid"
              FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -fbacktrace -finit-real=nan -Wimplicit-interface -Wtabs -Wall -Wextra"
              ;;
            enhanced)
              enable_debug_dbg="enhanced"
              FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -fbacktrace -finit-real=nan -Wimplicit-interface -Wtabs"
              ;;
            basic)
              enable_debug_dbg="basic"
              
              ;;
          esac   # [case: enable_debug, indent: 4, item: True]
          ;;
      esac   # [case: abi_fc_version, indent: 2, item: True]
      ;;
    open64)
      abi_fc_vendor_dbg="open64"
      abi_fc_version_dbg="default"
      abi_cpu_spec_dbg="default"
      case "${enable_debug}" in
        naughty)
          enable_debug_dbg="naughty"
          FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -g3 -fullwarn -trapuv -Wall -Wendif-labels -Wunused -ansi -fno-permissive -Wfloat-equal -Wunreachable-code -C"
          ;;
        verbose)
          enable_debug_dbg="verbose"
          
          ;;
        paranoid)
          enable_debug_dbg="paranoid"
          FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -g3 -fullwarn -trapuv -Wall -Wendif-labels -Wunused -ansi -fno-permissive -Wfloat-equal -Wunreachable-code"
          ;;
        enhanced)
          enable_debug_dbg="enhanced"
          FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -g3 -fullwarn -trapuv -Wall -Wendif-labels -Wunused"
          ;;
        basic)
          enable_debug_dbg="basic"
          
          ;;
      esac   # [case: enable_debug, indent: 2, item: True]
      ;;
    g95)
      abi_fc_vendor_dbg="g95"
      abi_fc_version_dbg="default"
      abi_cpu_spec_dbg="default"
      case "${enable_debug}" in
        naughty)
          enable_debug_dbg="naughty"
          FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -On -Sw -Ds"
          ;;
        verbose)
          enable_debug_dbg="verbose"
          
          ;;
        paranoid)
          enable_debug_dbg="paranoid"
          FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -On -Sw"
          ;;
        enhanced)
          enable_debug_dbg="enhanced"
          FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -On"
          ;;
        basic)
          enable_debug_dbg="basic"
          
          ;;
      esac   # [case: enable_debug, indent: 2, item: True]
      ;;
    fujitsu)
      abi_fc_vendor_dbg="fujitsu"
      abi_fc_version_dbg="default"
      abi_cpu_spec_dbg="default"
      case "${enable_debug}" in
        naughty)
          enable_debug_dbg="naughty"
          FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -AT -Du -Si"
          ;;
        verbose)
          enable_debug_dbg="verbose"
          
          ;;
        paranoid)
          enable_debug_dbg="paranoid"
          FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -AT -Du -Si"
          ;;
        enhanced)
          enable_debug_dbg="enhanced"
          FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -AT -Du -Si"
          ;;
        basic)
          enable_debug_dbg="basic"
          
          ;;
      esac   # [case: enable_debug, indent: 2, item: True]
      ;;
    pathscale)
      abi_fc_vendor_dbg="pathscale"
      abi_fc_version_dbg="default"
      abi_cpu_spec_dbg="default"
      case "${enable_debug}" in
        naughty)
          enable_debug_dbg="naughty"
          FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -trapuv -fullwarn -Wall -fno-permissive -Wendif-labels -ffortran-bounds-check"
          ;;
        verbose)
          enable_debug_dbg="verbose"
          
          ;;
        paranoid)
          enable_debug_dbg="paranoid"
          FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -trapuv -fullwarn -Wall -fno-permissive -Wendif-labels"
          ;;
        enhanced)
          enable_debug_dbg="enhanced"
          FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -trapuv -fullwarn -Wall"
          ;;
        basic)
          enable_debug_dbg="basic"
          
          ;;
      esac   # [case: enable_debug, indent: 2, item: True]
      ;;
    nag)
      abi_fc_vendor_dbg="nag"
      abi_fc_version_dbg="default"
      abi_cpu_spec_dbg="default"
      case "${enable_debug}" in
        naughty)
          enable_debug_dbg="naughty"
          FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -gline -mtrace=verbose -nan -info -C"
          FC_LDFLAGS_DEBUG="${FC_LDFLAGS_DEBUG} -mtrace=verbose"
          ;;
        verbose)
          enable_debug_dbg="verbose"
          
          ;;
        paranoid)
          enable_debug_dbg="paranoid"
          FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -gline -mtrace=verbose -nan -info"
          FC_LDFLAGS_DEBUG="${FC_LDFLAGS_DEBUG} -mtrace=verbose"
          ;;
        enhanced)
          enable_debug_dbg="enhanced"
          FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -gline -mtrace=verbose -nan"
          FC_LDFLAGS_DEBUG="${FC_LDFLAGS_DEBUG} -mtrace=verbose"
          ;;
        basic)
          enable_debug_dbg="basic"
          
          ;;
      esac   # [case: enable_debug, indent: 2, item: True]
      ;;
    intel)
      abi_fc_vendor_dbg="intel"
      abi_fc_version_dbg="default"
      abi_cpu_spec_dbg="default"
      case "${enable_debug}" in
        naughty)
          enable_debug_dbg="naughty"
          FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -debug all -check uninit -ftrapuv -traceback -warn all -fp-stack-check -check bounds -WB"
          ;;
        verbose)
          enable_debug_dbg="verbose"
          
          ;;
        paranoid)
          enable_debug_dbg="paranoid"
          FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -debug all -check uninit -ftrapuv -traceback -warn all -fp-stack-check"
          ;;
        enhanced)
          enable_debug_dbg="enhanced"
          FCFLAGS_DEBUG="${FCFLAGS_DEBUG} -debug all -check uninit -ftrapuv -traceback"
          ;;
        basic)
          enable_debug_dbg="basic"
          
          ;;
      esac   # [case: enable_debug, indent: 2, item: True]
      ;;
  esac   # [case: abi_fc_vendor, indent: 0, item: True]

  dnl Display settings
  AC_MSG_RESULT([${abi_fc_vendor_dbg}/${abi_fc_version_dbg}/${abi_cpu_spec_dbg}])

]) #ABI_FC_DBGFLAGS
