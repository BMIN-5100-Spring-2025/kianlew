
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
  bucket      = aws_s3_bucket.project_data_bucket.id
  acl         = "private"
}

resource "aws_s3_bucket_lifecycle_configuration" "project_data_bucket_expiration" {
  bucket = aws_s3_bucket.project_data_bucket.id

  rule {
    id     = "compliance-retention-policy"
    status = "Enabled"

    expiration {
      days = 100
    }
  }
}


resource "aws_ecr_repository" "lew_kianapp" {
  name = "lew_kianapp"

  image_scanning_configuration {
    scan_on_push = true
  }

  image_tag_mutability = "MUTABLE"

  tags = {
    Project = "lew_kianapp"
    Owner   = element(split("/", data.aws_caller_identity.current.arn), 1)
  }
}

resource "aws_cloudwatch_log_group" "lew_kianapp" {
  name              = "/ecs/lew_kianapp"
  retention_in_days = 30
}

import {
  to = aws_cloudwatch_log_group.lew_kianapp
  id = "/ecs/lew_kianapp"
}


resource "aws_ecs_task_definition" "lew_kianapp" {
  family                   = "lew_kianapp"
  requires_compatibilities = ["FARGATE"]
  cpu                      = 2048
  memory                   = 14336
  network_mode             = "awsvpc"
  execution_role_arn       = aws_iam_role.ecs_execution_role.arn
  task_role_arn            = aws_iam_role.ecs_task_role.arn

   ephemeral_storage {
     size_in_gib = 30
   }

  container_definitions = jsonencode([
    {
      name      = "lew_kianapp"
      image     = "061051226319.dkr.ecr.us-east-1.amazonaws.com/lew_kianapp:latest"
      cpu       = 2048
      memory    = 14336
      essential = true

      environment = [
        { name = "RUN_MODE",            value = "fargate" },
        { name = "S3_BUCKET_NAME",      value = "bmin5100-kianlew" },
        { name = "AWS_DEFAULT_REGION",  value = "us-east-1" },
        { name = "INPUT_DIR",           value = "/tmp/input" },
        { name = "OUTPUT_DIR",          value = "/tmp/output" }
      ]
    }
  ])
}


resource "aws_iam_role" "ecs_execution_role" {
  name = "ecs-execution-role"

  assume_role_policy = data.aws_iam_policy_document.ecs_task_trust.json
}

data "aws_iam_policy_document" "ecs_task_trust" {
  statement {
    actions = ["sts:AssumeRole"]
    principals {
      type        = "Service"
      identifiers = ["ecs-tasks.amazonaws.com"]
    }
  }
}

resource "aws_iam_role_policy_attachment" "ecs_execution_role_ecs_policy" {
  role       = aws_iam_role.ecs_execution_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy"
}

resource "aws_iam_role" "ecs_task_role" {
  name = "ecs-task-role"

  assume_role_policy = data.aws_iam_policy_document.ecs_task_trust.json
}

data "aws_s3_bucket" "my_imported_bucket" {
  bucket = "bmin5100-kianlew"
}

data "aws_iam_policy_document" "s3_access_for_task" {
  statement {
    actions = [
      "s3:GetObject",
      "s3:ListBucket",
      "s3:PutObject",
      "s3:DeleteObject",
    ]
    resources = [
      data.aws_s3_bucket.my_imported_bucket.arn,
      "${data.aws_s3_bucket.my_imported_bucket.arn}/*",
    ]
  }
}

resource "aws_iam_policy" "task_s3_policy" {
  name        = "task-s3-policy"
  description = "IAM policy for ECS task to access S3 bucket"
  policy      = data.aws_iam_policy_document.s3_access_for_task.json
}

resource "aws_iam_role_policy_attachment" "ecs_task_role_s3_policy" {
  role       = aws_iam_role.ecs_task_role.name
  policy_arn = aws_iam_policy.task_s3_policy.arn
}
