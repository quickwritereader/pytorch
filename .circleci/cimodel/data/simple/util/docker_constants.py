AWS_DOCKER_HOST = "308535385114.dkr.ecr.us-east-1.amazonaws.com"

# ARE YOU EDITING THIS NUMBER?  MAKE SURE YOU READ THE GUIDANCE AT THE
# TOP OF .circleci/config.yml
DOCKER_IMAGE_TAG = "8bdba785b1eac4d297d5f5930f979518012a56e0"


def gen_docker_image_path(container_type, container_tag=DOCKER_IMAGE_TAG):
    return "/".join([
        AWS_DOCKER_HOST,
        "pytorch",
        container_type + ":" + container_tag,
    ])


DOCKER_IMAGE_BASIC = gen_docker_image_path("pytorch-linux-xenial-py3.6-gcc5.4")

DOCKER_IMAGE_CUDA_10_2 = gen_docker_image_path("pytorch-linux-xenial-cuda10.2-cudnn7-py3-gcc7")

DOCKER_IMAGE_GCC7 = gen_docker_image_path("pytorch-linux-xenial-py3.6-gcc7")

def gen_mobile_docker_name(specifier):
    container_type = "pytorch-linux-xenial-py3-clang5-" + specifier
    return gen_docker_image_path(container_type)


DOCKER_IMAGE_ASAN = gen_mobile_docker_name("asan")

DOCKER_IMAGE_NDK = gen_mobile_docker_name("android-ndk-r19c")
