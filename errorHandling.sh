function catch_exit_status () {

    IFS=' ' read -a cmd <<< "$BASH_COMMAND"

    local exit_code=$1
    local lineNum=$2
    local src=$3

    set -- "${cmd[@]}"

    if [[ "${exit_code}" -eq 0 ]]
    then
        exit 0;
    fi

    if  [[ "${cmd[0]}" != "exit" ]]
    then
        echo -e >&2 "$(timestamp) ERROR: Error occurred in ${cmd[@]} in line $lineNum of $src (Exit code: $exit_code)."
    else
        echo -e >&2 "$(timestamp) ERROR: Error occurred in line $lineNum of $src (Exit code: $exit_code)."
    fi

    if [[ "${exit_code}" -gt 130 ]]
    then
        echo -e >&2 "$(timestamp) ERROR: Execution halted with a fatal error ($((exit_code - 128)))."
        exit ${exit_code}
    fi

    case ${exit_code} in
        2)
        echo -e >&2 "$(timestamp) ERROR: Missing argument."
        ;;
        3)
        echo -e >&2 "$(timestamp) ERROR: Missing input file."
        ;;
        4)
        echo -e >&2 "$(timestamp) ERROR: Invalid input file."
        ;;
        5)
        echo -e >&2 "$(timestamp) ERROR: Invalid BAM file."
        ;;
        6)
        echo -e >&2 "$(timestamp) ERROR: Invalid VCF file."
        ;;
        126)
        echo -e >&2 "$(timestamp) ERROR: Permission denied."
        ;;
        127)
        echo -e >&2 "$(timestamp) ERROR: Command not found. Pls check the prerequisites."
        ;;
        128)
        echo -e >&2 "$(timestamp) ERROR: Illegal arguments encountered."
        ;;
        129)
        echo -e >&2 "$(timestamp) ERROR: Execution halted with a fatal error (1)."
        ;;
        130)
        echo -e >&2 "$(timestamp) ERROR: Keyboard interrupt."
        ;;
        *)
        echo -e >&2 "$(timestamp) ERROR: Unknown error encountered."
        ;;
    esac

    exit ${exit_code}
}
