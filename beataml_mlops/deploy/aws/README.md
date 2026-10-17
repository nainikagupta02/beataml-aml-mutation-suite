# AWS deployment (not applied)

This Terraform module defines, but does not create, an AWS App Runner
service for `../../api`: an ECR repository, an IAM role App Runner uses to
pull from it, and the App Runner service itself (smallest instance size,
auto-deploy on new image push, `/health` as the health check path).

Nothing here has been run against any AWS account. Running it will create
real, billed resources. Review `main.tf` and the App Runner/ECR pricing page
before applying.

## Why App Runner

The API is a single stateless container with modest CPU/memory needs (a
~5,000-feature logistic-regression pipeline, no GPU). App Runner avoids
standing up a VPC, load balancer, and ECS cluster for that; if this ever
needs to share infrastructure with other services, moving to ECS Fargate
behind an ALB is a straightforward follow-up, not a rewrite of the image.

## Manual steps to actually deploy

1. Build and tag the image (see `../../api/README.md` for `prepare_artifacts.py` first):
   ```bash
   cd ../../api
   docker build -t beataml-mutation-api:latest .
   ```
2. `terraform init` and `terraform apply` once, just to create the ECR repo
   (the App Runner service will fail to start until an image exists at the
   tag it expects — see step 4 for the usual chicken-and-egg fix):
   ```bash
   cd ../deploy/aws
   terraform init
   terraform apply -target=aws_ecr_repository.beataml_api
   ```
3. Authenticate Docker to that ECR repo and push:
   ```bash
   aws ecr get-login-password --region <aws_region> \
     | docker login --username AWS --password-stdin <account_id>.dkr.ecr.<aws_region>.amazonaws.com
   docker tag beataml-mutation-api:latest <ecr_repository_url>:latest
   docker push <ecr_repository_url>:latest
   ```
4. Apply the rest:
   ```bash
   terraform apply
   ```
5. `terraform output service_url` gives the public HTTPS endpoint. Test it:
   ```bash
   curl https://<service_url>/health
   ```

## Teardown

```bash
terraform destroy
```

## Limitations

- No authentication/authorization in front of the API — App Runner's default
  URL is public. Add an API key check in `api/main.py` or put API Gateway /
  a WAF in front before exposing this beyond local testing.
- No autoscaling, logging export, or alarm configuration beyond App Runner's
  defaults.
- This is a feasibility-study model (see `../../README.md` and the parent
  projects' READMEs); nothing about deploying it changes that.
