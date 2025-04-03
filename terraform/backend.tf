terraform {
  backend "s3" {
    bucket         = "bmin5100-terraform-state"
    key            = "kian.lew@pennmedicine.upenn.edu-template/terraform.tfstate"
    region         = "us-east-1"
    encrypt        = true
  }
}