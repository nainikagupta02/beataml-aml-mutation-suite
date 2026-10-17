# AWS App Runner deployment for the BEAT-AML mutation-probability API.
#
# NOT APPLIED. This defines the infrastructure; it does not run `terraform
# apply` against any account. See README.md in this directory for the manual
# steps (build/push the image, then plan/apply) and for cost before you do.

terraform {
  required_version = ">= 1.5"
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 5.0"
    }
  }
}

provider "aws" {
  region = var.aws_region
}

resource "aws_ecr_repository" "beataml_api" {
  name                 = var.ecr_repository_name
  image_tag_mutability = "MUTABLE"
  image_scanning_configuration {
    scan_on_push = true
  }
}

resource "aws_ecr_lifecycle_policy" "beataml_api" {
  repository = aws_ecr_repository.beataml_api.name
  policy = jsonencode({
    rules = [{
      rulePriority = 1
      description  = "Keep last 10 images"
      selection = {
        tagStatus   = "any"
        countType   = "imageCountMoreThan"
        countNumber = 10
      }
      action = { type = "expire" }
    }]
  })
}

# App Runner needs an IAM role with permission to pull from this private ECR repo.
resource "aws_iam_role" "apprunner_ecr_access" {
  name = "${var.service_name}-ecr-access"
  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect    = "Allow"
      Principal = { Service = "build.apprunner.amazonaws.com" }
      Action    = "sts:AssumeRole"
    }]
  })
}

resource "aws_iam_role_policy_attachment" "apprunner_ecr_access" {
  role       = aws_iam_role.apprunner_ecr_access.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AWSAppRunnerServicePolicyForECRAccess"
}

resource "aws_apprunner_service" "beataml_api" {
  service_name = var.service_name

  source_configuration {
    authentication_configuration {
      access_role_arn = aws_iam_role.apprunner_ecr_access.arn
    }
    image_repository {
      image_identifier      = "${aws_ecr_repository.beataml_api.repository_url}:${var.image_tag}"
      image_repository_type = "ECR"
      image_configuration {
        port = "8000"
      }
    }
    # Re-deploy automatically when a new image is pushed with the same tag.
    auto_deployments_enabled = true
  }

  instance_configuration {
    cpu    = var.instance_cpu    # "0.25 vCPU" smallest tier; the model is small, no GPU needed
    memory = var.instance_memory # "0.5 GB" is enough for a ~5,000-feature logistic-regression pipeline
  }

  health_check_configuration {
    protocol            = "HTTP"
    path                = "/health"
    interval            = 10
    timeout             = 5
    healthy_threshold   = 1
    unhealthy_threshold = 3
  }

  tags = {
    Project = "beataml-mlops"
    Purpose = "frozen-model-serving"
  }
}
