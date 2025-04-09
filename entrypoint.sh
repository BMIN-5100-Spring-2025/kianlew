#!/bin/bash
echo "Container starting with ENVIRONMENT=$ENVIRONMENT, S3_BUCKET=$S3_BUCKET, AWS_DEFAULT_REGION=$AWS_DEFAULT_REGION"

exec Rscript /usr/src/app/app/final_project_template_bmin5100.R