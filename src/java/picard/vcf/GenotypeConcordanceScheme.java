package picard.vcf;

import htsjdk.samtools.util.StringUtil;
import picard.PicardException;
import picard.vcf.GenotypeConcordanceStates.*;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * This defines for each valid TruthState and CallState tuple, the set of contingency table entries that to which the tuple should contribute.
 * @author nhomer
 */
public abstract class GenotypeConcordanceScheme{

    /** The underlying scheme */
    protected final Map<TruthAndCallStates, ContingencyState[]> scheme = new HashMap<TruthAndCallStates, ContingencyState[]>();

    /** These are convenience variables for defining a scheme.  NA means that such a tuple should never be observed. */
    public static final ContingencyState[]    NA       = {ContingencyState.NA};
    protected static final ContingencyState[] EMPTY    = {ContingencyState.EMPTY};
    protected static final ContingencyState[] TP_ONLY  = {ContingencyState.TP};
    protected static final ContingencyState[] FP_ONLY  = {ContingencyState.FP};
    protected static final ContingencyState[] TN_ONLY  = {ContingencyState.TN};
    protected static final ContingencyState[] FN_ONLY  = {ContingencyState.FN};
    protected static final ContingencyState[] TP_FN    = {ContingencyState.TP, ContingencyState.FN};
    protected static final ContingencyState[] TP_FP    = {ContingencyState.TP, ContingencyState.FP};
    protected static final ContingencyState[] TP_TN    = {ContingencyState.TP, ContingencyState.TN};
    protected static final ContingencyState[] FP_FN    = {ContingencyState.FP, ContingencyState.FN};
    protected static final ContingencyState[] FP_TN    = {ContingencyState.FP, ContingencyState.TN};
    protected static final ContingencyState[] FP_TN_FN = {ContingencyState.FP, ContingencyState.TN, ContingencyState.FN};
    protected static final ContingencyState[] TP_FP_FN = {ContingencyState.TP, ContingencyState.FP, ContingencyState.FN};
    protected static final ContingencyState[] TN_FN    = {ContingencyState.TN, ContingencyState.FN};

    /** Has this scheme been previously validated */
    private boolean isValidated = false;

    /** Initiates the scheme by adding rows to it with calls based on the scheme design
     *  Current schemes include
     *  GA4GHScheme
     *  GA4GHSchemeWithMissingAsNoCall
     *  */
    public abstract void initiateScheme();

    /**
     * Adds a row to the scheme
     * @param callState the call state (row)
     * @param concordanceStateArrays the concordance state arrays for each truth value, in order
     */
    protected void addRow(final CallState callState, final ContingencyState[]... concordanceStateArrays) {
        if (concordanceStateArrays.length != TruthState.values().length) {
            throw new PicardException("Length mismatch between concordanceStateArrays and TruthState.values()");
        }
        for (int i = 0; i < concordanceStateArrays.length; i++) {
            scheme.put(new TruthAndCallStates(TruthState.values()[i], callState), concordanceStateArrays[i]);
        }
    }

    /**
     * Get the concordance state array associate with the given truth state and call state tuple.
     */
    public ContingencyState[] getConcordanceStateArray(final TruthState truthState, final CallState callState) {
        return this.getConcordanceStateArray(new TruthAndCallStates(truthState, callState));
    }

    /**
     * Get the concordance state array associate with the given truth state and call state tuple.
     */
    public ContingencyState[] getConcordanceStateArray(final TruthAndCallStates truthAndCallStates) {
        return this.scheme.get(truthAndCallStates);
    }

    /**
     * Get the contingency state array as a parse-able string
     */
    public String getContingencyStateString(final TruthState truthState, final CallState callState) {
        final ContingencyState[] contingencyStateArray = getConcordanceStateArray(truthState, callState);
        return (contingencyStateArray.length == 0) ? "EMPTY" : StringUtil.join(",", contingencyStateArray);
    }

    /**
     * Get the contingency state array as a set
     * @param contingencyStateArray
     * @return
     */
    public Set<ContingencyState> getContingencyStateSet(final ContingencyState[] contingencyStateArray) {
        final Set<ContingencyState> contingencyStateSet = new HashSet<ContingencyState>();
        Collections.addAll(contingencyStateSet, contingencyStateArray);
        return contingencyStateSet;
    }


    /**
     * Check that all cells in the scheme exist.
     * @throws PicardException if a missing tuple was found.
     */
    public void validateScheme() throws PicardException {
        if (!isValidated) {
            for (final TruthState truthState : TruthState.values()) {
                for (final CallState callState : CallState.values()) {
                    if (!scheme.containsKey(new TruthAndCallStates(truthState, callState))) {
                        throw new PicardException(String.format("Missing scheme tuple: [%s, %s]", truthState.name(), callState.name()));
                    }
                }
            }
        }

        isValidated = true;
    }
}

