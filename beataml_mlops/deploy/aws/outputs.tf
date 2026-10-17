output "service_url" {
  description = "Public HTTPS URL App Runner assigns once the service is running."
  value       = aws_apprunner_service.beataml_api.service_url
}

output "ecr_repository_url" {
  description = "Push images here: docker push <ecr_repository_url>:<tag>"
  value       = aws_ecr_repository.beataml_api.repository_url
}
