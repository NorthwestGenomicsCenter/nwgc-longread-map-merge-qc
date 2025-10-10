public class Utils {

    public static defineReadGroup(fastq, sampleId) {

        // Create the tags
        def readGroupTags = []
        readGroupTags.add("ID:" + fastq)
        readGroupTags.add("PL:PACBIO")
        readGroupTags.add("PU:" + fastq)
        readGroupTags.add("SM:" + sampleId)

        // Create the read group to the
        return "'@RG\\t" + readGroupTags.join("\\t") + "'"
    }

}
