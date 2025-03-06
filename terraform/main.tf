resource "aws_s3_bucket" "project_data_bucket" {
  bucket = "bmin5100-kianlew"

  tags = {
    Owner = element(split("/", data.aws_caller_identity.current.arn), 1)
  }
}

resource "aws_s3_bucket_ownership_controls" "project_data_bucket_ownership_controls" {
  bucket = aws_s3_bucket.project_data_bucket.id
  rule {
    object_ownership = "BucketOwnerPreferred"
  }
}

resource "aws_s3_bucket_acl" "project_data_bucket_acl" {
  depends_on = [aws_s3_bucket_ownership_controls.project_data_bucket_ownership_controls]

  bucket = aws_s3_bucket.project_data_bucket.id
  acl    = "private"
}

resource "aws_s3_bucket_lifecycle_configuration" "project_data_bucket_expiration" {
  bucket = aws_s3_bucket.project_data_bucket.id

  rule {
    id      = "compliance-retention-policy"
    status  = "Enabled"

    expiration {
	  days = 100
    }
  }
}