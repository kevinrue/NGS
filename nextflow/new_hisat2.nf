

params.weblink = "ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.5-Linux_x86_64.zip"

Channel
	.value(params.weblink)
	.set { weblinkChannel }

process downloadLinuxBinary {
    input:
    val weblinkChannel
    
    output:
    file "archive" into archive
      
    script:
    //
    // Download the new 
    //
    """
    echo wget weblinkChannel
    """
}

/*
 * print the channel content
 */
archive.subscribe { println it }
