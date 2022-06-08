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
    
    while :
    do
        read -p "$(timestamp) PROMPT: Do you want to report the error via email? (y/n)" -n 1 send_email
        if [ "${send_email}" == "y" ]
        then
            echo -e "$(timestamp) INFO: An email will be sent to report the incident."
            emailAddr='u3005579@connect.hku.hk'
            temp_dir="/home/louisshe/ngs_scripts_louis/tmp"
            echo -e "Dear Yangyxt,\n" > "${temp_dir}"/error_$(date +%Y%m%d).txt
            echo -e "An error is reported on $(date +%F\ %T) with exit code ${exit_code}.\n" >> "${temp_dir}"/error_$(date +%Y%m%d).txt
            echo -e "Below is the list of log files found in ~/tmp/ for your reference.\n" >> "${temp_dir}"/error_$(date +%Y%m%d).txt
            find ${temp_dir}/*.log -type f >> "${temp_dir}"/error_$(date +%Y%m%d).txt
            echo -e '\n'
            mail -s 'Error running pipelines $(date +%Y%m%d)' ${emailAddr} < "${temp_dir}"/error_$(date +%Y%m%d).txt || { echo -e >&2 "$(timestamp) ERROR: Unable to send email!" && exit 129; }
            break;
        elif [ "${send_email}" == "n" ]
        then
            echo -e '\n'
            break;
        fi
        echo -e '\n'
    done
    exit ${exit_code}
}
