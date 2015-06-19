package picard.vcf;

/**
 * Created by kbergin on 6/19/15.
 */
public class GenotypeConcordanceSchemeFactory {
    public GenotypeConcordanceScheme getScheme(final boolean isMissingHomRef){
        if(isMissingHomRef){
            return new GA4GHSchemeWithMissingAsNoCall();
        }
        else{
            return new GA4GHScheme();
        }
    }
}
