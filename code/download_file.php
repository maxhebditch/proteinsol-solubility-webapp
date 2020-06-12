<?php

$app =  $_GET["app"];
$timestamp =  $_GET["timestamp"];
$file =  $_GET["file"];
$idname = $_GET["idname"];
$dirname = $_GET["dirname"];

$filename = "/var/www/html/results/run-" . $timestamp . "/" . $file;

$cleanname = $idname . "_" . $file;

header("Pragma: public");
header("Expires: 0");
header("Cache-Control: must-revalidate, post-check=0, pre-check=0");
header("Cache-Control: public");
header("Content-Description: File Transfer");
header("Content-type: application/octet-stream");
header("Content-Disposition: attachment; filename=\"".$cleanname."\"");
header("Content-Transfer-Encoding: binary");
header("Content-Length: ".filesize($filename));
ob_end_flush();
@readfile($filename);

?>
