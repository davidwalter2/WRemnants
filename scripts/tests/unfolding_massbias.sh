# unfolding bias test
# 0) Generate card with nominal model and pseudodata sets with different masses 
# 1) Unfold pseudodata with mW set to values of {0, 10, 20, ...}
# 2) Generate card with nominal theory model
# 3) Fit theory model to unfolded pseudo data 
# 4) Compare observed mass pulls with pseudodata value

# hyperparameters
# a) uncertainty on mass when unfolding
# a1) None
# a2) 100MeV
# b) Histmaker type
# b1) mz_wlike
# b2) mw
# c) transverse mass cut (only for mz_wlike possible on real data)
# c1) with cuts
# c2) without cut 

# settings
LABEL=Z
ANALYSIS=ZMassWLike
STATONLY=_statOnly
OPT="--ewUnc "
if [ -n "$STATONLY" ]; then
    OPT="--doStatOnly"
fi

CMSSW_BASE=/home/${USER:0:1}/${USER}/CMSSW_10_6_30/src/
HISTMAKER_FILE=/scratch/${USER}/results_histmaker/231106_unfolding_wlike/mz_wlike_with_mu_eta_pt_scetlib_dyturboCorr_unfolding_mtCut0.hdf5
COMBINE_OUTDIR=/scratch/${USER}/CombineStudies/unfolding_massbias/231109_wlike_mtCut0_mZUnc100MeV${STATONLY}

COMBINE_ANALYSIS_OUTDIR=${COMBINE_OUTDIR}/${ANALYSIS}_eta_pt_charge${STATONLY}/
COMBINE_ANALYSIS_PATH=${COMBINE_ANALYSIS_OUTDIR}/${ANALYSIS}.hdf5

# 0) Generate card with nominal model and pseudodata sets with different masses 
if [ -e $COMBINE_ANALYSIS_PATH ]; then
    echo "The file $COMBINE_ANALYSIS_PATH exists, continue using it."
else
    echo "The file does not exists, produce it."
    ./scripts/ci/run_with_singularity.sh scripts/ci/setup_and_run_python.sh scripts/combine/setupCombine.py \
        -i $HISTMAKER_FILE -o $COMBINE_OUTDIR --hdf5 --sparse --unfolding \
        --pseudoData massWeight$LABEL --pseudoDataAxes massShift --pseudoDataIdxs -1 --ewUnc $OPT
fi

# 1) Unfold pseudodata with mW set to values of {0, 10, 20, ...}
for ((i=0; i<=20; i+=5)); do
    IABS=$(( (i-10) * 10 ))
    if [ $IABS -lt 0 ]; then
        UPDOWN="Down"
        IABS=$((-1*$IABS))
    elif [ $IABS -eq 0 ]; then
        UPDOWN=""
    else
        UPDOWN="Up"
    fi
    BIN=massShift${LABEL}${IABS}MeV$UPDOWN
    PSEUDO="massWeight${LABEL}_massShift_${BIN}"

    echo "Perform unfolding with index = $i : $PSEUDO"

    # 1) Unfold pseudodata
    FITRESULT=${COMBINE_ANALYSIS_OUTDIR}/fitresults_123456789_${BIN}.hdf5

    if [ -e $FITRESULT ]; then
        echo "The file $FITRESULT exists, continue using it."
    else
        cmssw-cc7 --command-to-run scripts/ci/setup_and_run_combine.sh $CMSSW_BASE $COMBINE_ANALYSIS_OUTDIR \
            ${ANALYSIS}.hdf5 -p $PSEUDO --postfix $BIN --binByBinStat --correlateXsecStat #--doImpacts 
    fi

    # 2)  Generate card with nominal theory model
    THEOMODEL_DIR=${COMBINE_OUTDIR}/${ANALYSIS}_qGen_ptGen_absEtaGen${STATONLY}_${BIN}/
    THEOMODEL=${ANALYSIS}.hdf5
    THEOMODEL_PATH=${THEOMODEL_DIR}/${THEOMODEL}

    if [ -e $THEOMODEL_PATH ]; then
        echo "The file $THEOMODEL_PATH exists, continue using it."
    else
        ./scripts/ci/run_with_singularity.sh scripts/ci/setup_and_run_python.sh scripts/combine/setupCombine.py \
            -i $HISTMAKER_FILE -o $COMBINE_OUTDIR --hdf5 --fitvar qGen-ptGen-absEtaGen --postfix $BIN $OPT \
            --fitresult $FITRESULT
    fi

    # 3) Fit theory model to unfolded pseudo data 
    FITRESULT_FINAL=${THEOMODEL_DIR}/fitresults_123456789.hdf5

    if [ -e $FITRESULT_FINAL ]; then
        echo "The file $FITRESULT_FINAL exists, continue using it."
    else
        cmssw-cc7 --command-to-run scripts/ci/setup_and_run_combine.sh $CMSSW_BASE $THEOMODEL_DIR \
            $THEOMODEL --chisqFit --externalCovariance # --doImpacts 
    fi

    # 4) Compare observed mass pulls with pseudodata value


done




