<?php 

/* adds in generic header */
header('Content-type: text/html; charset=utf-8');
include 'header.php';

// Turn off output buffering
ini_set('output_buffering', 'off');
         
//Flush (send) the output buffer and turn off output buffering
while (@ob_end_flush());
         
// Implicitly flush the buffer(s)
ini_set('implicit_flush', true);
ob_implicit_flush(true);

#function to run the loading wheel
function showint($timeStamp) {
	ob_start();
	?>
	<div id="INTloading">
	<div style="text-align: center;height:100vh">
	<br><br>
	<h3>Loading...</h3>
	<br><br>
	<p>Calculation: solubility</p>
	<p>Job id = <?php echo $timeStamp ?></p>
	<br><br>
	<div class="loader"></div>
	</div>
	</div>
	<?php 
	ob_flush();
	flush();
}

#if sequence insertion is empty error
if (empty($_POST['sequence-input'])) {
        echo "<p>It appears you didn't input any sequence. Please input an amino acid sequence for calculation</p>";
        echo "<br>";
        echo "<p>For example</p>";
        echo "<div style='word-wrap:break-word;overflow:hidden'>";
        echo "<code>&gt P00547</code>";
        echo "<br>";
        echo "<code>MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGRFADKLPSEPRENIVYQCWERFCQELGKQIPVAMTLEKNMPIGSGLGSSACSVVAALMAMNEHCGKPLNDTRLLALMGELEGRISGSIHYDNVAPCFLGGMQLMIEENDIISQQVPGFDEWLWVLAYPGIKVSTAEARAILPAQYRRQDCIAHGRHLAGFIHACYSRQPELAAKLMKDVIAEPYRERLLPGFRQARQAVAEIGAVASGISGSGPTLFALCDKPETAQRVADWLGKNYLQNQEGFVHICRLDTAGARVLEN</code>";
        echo "<br>";
        echo "</div>";
        echo "<br>";
        echo "<p>Please contact us at protein-sol@manchester.ac.uk if you think it should have worked.</p>";
        exit();
}

/* pulls in the sequence and converts to upper case and removes whitespace */
$sequenceInput = $_POST['sequence-input'];

#assign random timestamp
require_once "random.php";
$timeStamp = bin2hex(random_bytes(10));
#assign dir and filename
$tempdir = "/tmp/run-" . $timeStamp . "/";

#make temporary dir
mkdir($tempdir);

#base filename
$basefileName = "sequenceIn-" . $timeStamp;
#filename to save uploaded sequence
$fileName = $tempdir . "sequenceIn-" . $timeStamp;

#save sequence to file
$sequenceSave = fopen($fileName,"w");
fwrite($sequenceSave,$sequenceInput);


#code to tidy the fasta input that will generate an error file if there were any problems
#sequences that fail to run will be stuck in /tmp/run-TIMESTAMP until they are cleaned out
#by the cronjob
shell_exec("python3 /var/www/html/code/validate_fasta_write_file.py --input $fileName --min 21 --max 10000");
$errorfileName = $fileName . "-ERROR";

#if the error file exists, there must have been some errors in the submitted sequence
if (file_exists($errorfileName)){
        #open the list of errors as an array
	$error_array = file($errorfileName);
	echo '<div style="text-align: left;">';
	echo "JOB ID: $timeStamp";
	echo "<br>";
	echo "<br>";
	echo "<p>The following error(s) was identified in the submitted sequence:</p><br>";
	echo "<div style='padding-left:1em;padding-right:1em;'>";
        #iterate through the list of errors and write them to html
	foreach($error_array as $line){
			echo "<p id='codelike'>$line</p>";
	}
	echo "</div>";
	echo "<br>";
	echo "</div>";
        exit();
}

#END ERROR CHECKING INPUT FILE

#Get filename
$idname = shell_exec("grep '>' $fileName");
$idname = ltrim($idname,"> ");
//remove eol
$idname = trim(preg_replace('/\s+/', ' ', $idname));

#run the loading wheel
showint($timeStamp);
ob_flush();
flush();

#move to the temporary directory
#run jim's script
shell_exec("/var/www/html/code/run.sh $fileName $timeStamp '$idname'");

#if the prediction file was made correctly
if (file_exists("/tmp/run-$timeStamp/seq_prediction.txt")) {
                #copy files to the results dir
                shell_exec("mv '/tmp/run-$timeStamp' '/var/www/html/results/'");
                
?>
		<script type="text/javascript">
		    var timestamp = "<?php echo $timeStamp ?>";
		    window.location = "/results/run-"+timestamp+"/results.html";
		</script>

<?php
}else{
		echo"<style type='text/css'>#INTloading{";
		echo"display:none;";
		echo"}</style>";
                echo "<p>Something unexpected went wrong during calculation</p>";
                echo "<br>";
                echo "<p>Please input a single amino acid sequence for calculation";
                echo "<br>";
                echo "<p>For example</p>";
                echo "<div style='word-wrap:break-word;overflow:hidden'>";
                echo "<code>&gt P00547</code>";
                echo "<br>";
                echo "<code>MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGRFADKLPSEPRENIVYQCWERFCQELGKQIPVAMTLEKNMPIGSGLGSSACSVVAALMAMNEHCGKPLNDTRLLALMGELEGRISGSIHYDNVAPCFLGGMQLMIEENDIISQQVPGFDEWLWVLAYPGIKVSTAEARAILPAQYRRQDCIAHGRHLAGFIHACYSRQPELAAKLMKDVIAEPYRERLLPGFRQARQAVAEIGAVASGISGSGPTLFALCDKPETAQRVADWLGKNYLQNQEGFVHICRLDTAGARVLEN</code>";
                echo "<br>";
                echo "</div>";
                echo "<br>";
                echo "<p>Please contact us at protein-sol@manchester.ac.uk if you think it should have worked.</p>";
                exit();
}




?>
