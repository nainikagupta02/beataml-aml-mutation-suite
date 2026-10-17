variable "aws_region" {
  description = "AWS region to deploy into."
  type        = string
  default     = "us-east-1"
}

variable "service_name" {
  description = "App Runner service name."
  type        = string
  default     = "beataml-mutation-api"
}

variable "ecr_repository_name" {
  description = "ECR repository name for the API image."
  type        = string
  default     = "beataml-mutation-api"
}

variable "image_tag" {
  description = "Image tag to deploy. Push this tag to ECR before applying."
  type        = string
  default     = "latest"
}

variable "instance_cpu" {
  description = "App Runner instance CPU size."
  type        = string
  default     = "0.25 vCPU"
}

variable "instance_memory" {
  description = "App Runner instance memory size."
  type        = string
  default     = "0.5 GB"
}
